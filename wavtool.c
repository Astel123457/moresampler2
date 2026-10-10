#include "wavtool.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ciglet/ciglet.h>
#include <libllsm/llsm.h>

static char* build_debug_log_path(const char* filename) {
    size_t len = strlen(filename) + 14;
    char* path = (char*)malloc(len);
    if (!path) return NULL;
    snprintf(path, len, "%s", filename);
    char* ext = strrchr(path, '.');
    if (ext) strcpy(ext, ".ms2debug.log");
    else strcat(path, ".ms2debug.log");
    return path;
}

static char* build_session_path(const char* filename) {
    size_t len = strlen(filename) + 16;
    char* path = (char*)malloc(len);
    if (!path) return NULL;
    snprintf(path, len, "%s", filename);
    char* ext = strrchr(path, '.');
    if (ext) strcpy(ext, ".ms2session");
    else strcat(path, ".ms2session");
    return path;
}

static void append_debug_log(const char* path, int reset, const char* fmt, ...) {
    if (!path) return;
    FILE* f = fopen(path, reset ? "w" : "a");
    if (!f) return;
    va_list args;
    va_start(args, fmt);
    vfprintf(f, fmt, args);
    va_end(args);
    fclose(f);
}

static int can_resume_output_wav(const char* filename) {
    unsigned char header[12];
    FILE* f = fopen(filename, "rb");
    if (!f) return 0;
    size_t nread = fread(header, 1, sizeof(header), f);
    fclose(f);
    if (nread < sizeof(header)) return 0;
    return memcmp(header, "RIFF", 4) == 0 && memcmp(header + 8, "WAVE", 4) == 0;
}

static int session_contains_input(const char* session_path, const char* input_path) {
    if (!session_path || !input_path) return 0;
    FILE* f = fopen(session_path, "rb");
    if (!f) return 0;

    char line[4096];
    int found = 0;
    while (fgets(line, sizeof(line), f)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }
        if (strcmp(line, input_path) == 0) {
            found = 1;
            break;
        }
    }
    fclose(f);
    return found;
}

static void session_append_input(const char* session_path, const char* input_path) {
    if (!session_path || !input_path) return;
    FILE* f = fopen(session_path, "a");
    if (!f) return;
    fprintf(f, "%s\n", input_path);
    fclose(f);
}

static void reset_phrase_session(const char* output_llsm, const char* session_path) {
    if (output_llsm) remove(output_llsm);
    if (session_path) remove(session_path);
}

// Apply a volume scalar (0..100) to a single L1 frame's VTMAGN.
static void apply_volume_to_frame(llsm_container* frame, float vol) {
    if (vol <= 0.0f) {
        FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
        if (f0) f0[0] = 0.0f;
        return;
    }
    if (fabsf(vol - 100.0f) < 1e-6f) return;
    FP_TYPE* vt = (FP_TYPE*)llsm_container_get(frame, LLSM_FRAME_VTMAGN);
    if (!vt) return;
    int nspec = llsm_fparray_length(vt);
    float vol_db = 20.0f * log10f(vol / 100.0f);
    for (int j = 0; j < nspec; j++)
        vt[j] += vol_db;
}

// Build a per-frame envelope volume array (0..100) from the wavtool envelope points.
// Timeline: (0,0) → (p1,v1) → (p2,v2) [→ (p5,v5)] → (total-p3,v3) → (total-p4,v4)
static float* build_envelope(const wavtool_data* d, int nfrm, float ms_per_frame) {
    float* env = (float*)malloc(sizeof(float) * (size_t)nfrm);
    if (!env) return NULL;

    float total_ms = (float)nfrm * ms_per_frame;

    typedef struct { float t; float v; } EP;
    EP pts[7];
    int npts = 0;

    pts[npts++] = (EP){ 0.0f,                         d->v1 };
    pts[npts++] = (EP){ d->p1,                        d->v1 };
    // p2 is a duration after p1, not an absolute position
    pts[npts++] = (EP){ d->p1 + d->p2,               d->v2 };
    if (d->p5 > 0.0f)
        pts[npts++] = (EP){ d->p5,                   d->v5 };
    // p3 = duration of the v3→v4 fade; p4 = duration of the v4→0 final fade
    pts[npts++] = (EP){ total_ms - d->p4 - d->p3,    d->v3 };
    pts[npts++] = (EP){ total_ms - d->p4,             d->v4 };
    pts[npts++] = (EP){ total_ms,                     0.0f  };

    for (int i = 0; i < nfrm; i++) {
        float t = (float)i * ms_per_frame;
        float vol = 100.0f;
        for (int j = 0; j < npts - 1; j++) {
            if (t <= pts[j + 1].t) {
                float span = pts[j + 1].t - pts[j].t;
                float alpha = (span > 1e-6f) ? (t - pts[j].t) / span : 1.0f;
                if (alpha < 0.0f) alpha = 0.0f;
                if (alpha > 1.0f) alpha = 1.0f;
                vol = pts[j].v + alpha * (pts[j + 1].v - pts[j].v);
                break;
            }
        }
        env[i] = vol;
    }
    return env;
}

static void apply_gain_db_to_frame(llsm_container* frame, float gain_db) {
    FP_TYPE* vt = (FP_TYPE*)llsm_container_get(frame, LLSM_FRAME_VTMAGN);
    if (!vt) return;

    int nspec = llsm_fparray_length(vt);
    for (int j = 0; j < nspec; j++)
        vt[j] += gain_db;
}

int wavtool(wavtool_data* data) {
    // 1. Load input L1 chunk from the resampler's .llsm2 cache.
    //    No fallback to re-analysis — the resampler must have created this.
    char* input_llsm = build_llsm_path(data->input);
    if (!input_llsm) {
        fprintf(stderr, "wavtool: out of memory building input cache path\n");
        return 1;
    }
    int nfrm_in = 0, fs = 44100, nbit = 16;
    llsm_chunk* input = NULL;
    float alias_gain_db = 0.0f;
    float out_alias_gain_db = 0.0f;
    {
        FILE* f = fopen(input_llsm, "rb");
        if (f) { fclose(f); input = read_llsm(input_llsm, &nfrm_in, &fs, &nbit, &alias_gain_db); }
    }
    free(input_llsm);
    if (!input || nfrm_in == 0) {
        fprintf(stderr, "wavtool: no L1 cache found for input '%s' — "
                        "resampler must run before wavtool\n", data->input);
        if (input) llsm_delete_chunk(input);
        return 1;
    }

    // 2. Load output L1 chunk. Missing = first note for this render, not an error.
    //    If the output .wav doesn't exist, OpenUTAU started a fresh render — any
    //    stale .llsm2 cache from a previous render must be ignored.
    char* output_llsm = build_llsm_path(data->output);
    if (!output_llsm) {
        fprintf(stderr, "wavtool: out of memory building output cache path\n");
        llsm_delete_chunk(input);
        return 1;
    }
    char* session_path = build_session_path(data->output);
    int nfrm_out = 0, fs_out = 0, nbit_out = 0;
    llsm_chunk* out_chunk = NULL;
    {
        // Only resume the accumulated phrase cache if the output WAV already
        // exists as a readable WAV file. A fresh rerender may leave behind an
        // empty or invalid file path at the same location, which must not
        // revive a stale .llsm2 phrase cache from the previous render.
        int resume_output = data->use_position
            ? !data->first
            : can_resume_output_wav(data->output);
        int repeated_input = !data->use_position &&
            session_contains_input(session_path, data->input);
        if (!resume_output || repeated_input) {
            reset_phrase_session(output_llsm, session_path);
        } else if (resume_output) {
            FILE* f = fopen(output_llsm, "rb");
            if (f) { fclose(f); out_chunk = read_llsm(output_llsm, &nfrm_out, &fs_out, &nbit_out, &out_alias_gain_db); }
        }
    }
    if (out_chunk && fs_out > 0) { fs = fs_out; nbit = nbit_out; }

    char* debug_log = build_debug_log_path(data->output);

    // 3. ms → frames  (nhop=128 assumed, matching the resampler)
    float ms_per_frame = (float)(128 * 1000) / (float)fs;
    int offset_frames = (int)roundf(data->offset / ms_per_frame);
    int ovr_frames    = (int)roundf(data->ovr    / ms_per_frame);
    printf("[wavtool] input=%s\n", data->input);
    printf("[wavtool] output=%s\n", data->output);
    printf("[wavtool] length=%.3f ms  offset(stp)=%.3f ms  ovr=%.3f ms\n",
           data->length, data->offset, data->ovr);
    printf("[wavtool] ms_per_frame=%.4f  nfrm_in=%d  offset_frames=%d  ovr_frames=%d\n",
           ms_per_frame, nfrm_in, offset_frames, ovr_frames);
    if (ovr_frames < 0) ovr_frames = 0;

    int first_note = (nfrm_out == 0);
    float effective_offset_ms = data->offset;
    int effective_offset_frames = (int)roundf(effective_offset_ms / ms_per_frame);
    if (effective_offset_frames < 0) effective_offset_frames = 0;
    if (effective_offset_frames >= nfrm_in) effective_offset_frames = nfrm_in - 1;
    int delta_frames = (int)roundf(data->delta / ms_per_frame);

    int length_frames = (int)roundf(data->length / ms_per_frame);
    if (length_frames < 1) length_frames = 1;
    int in_start = effective_offset_frames;
    int in_count = nfrm_in - in_start;
    if (in_count <= 0) {
        fprintf(stderr, "wavtool: no input frames remain after effective offset trim\n");
        llsm_delete_chunk(input);
        if (out_chunk) llsm_delete_chunk(out_chunk);
        free(output_llsm);
        free(session_path);
        return 1;
    }
    int write_frames = length_frames;
    if (write_frames < 1) write_frames = 1;
    if (in_count > write_frames) in_count = write_frames;
    printf("[wavtool] in_start=%d  in_count=%d  length_frames=%d  ovr_frames=%d  nfrm_out=%d\n",
           in_start, in_count, length_frames, ovr_frames, nfrm_out);

    // 4. Per-frame envelope (volume 0..100)
    float* env = build_envelope(data, in_count, ms_per_frame);
    if (!env) {
        fprintf(stderr, "wavtool: out of memory for envelope\n");
        llsm_delete_chunk(input);
        if (out_chunk) llsm_delete_chunk(out_chunk);
        free(output_llsm);
        free(session_path);
        return 1;
    }

    int eff_ovr = (out_chunk && nfrm_out > 0) ? ovr_frames : 0;
    if (eff_ovr > nfrm_out) eff_ovr = nfrm_out;
    if (eff_ovr > in_count) eff_ovr = in_count;

    // 5. Concatenate: match sillyseam's effective overlap behavior by only
    // removing as many old frames as can actually be overlapped.
    int absolute_start = data->use_position
        ? (int)roundf(data->position / ms_per_frame) : 0;
    if (absolute_start < 0) absolute_start = 0;
    int keep_out = first_note ? 0 : (nfrm_out - eff_ovr);
    if (keep_out < 0) keep_out = 0;
    int new_nfrm = data->use_position
        ? max(nfrm_out, absolute_start + in_count)
        : keep_out + in_count;

    append_debug_log(debug_log, first_note,
        "[wavtool]\n"
        "  input=%s\n"
        "  output=%s\n"
        "  length_ms=%.3f offset_ms=%.3f delta_ms=%.3f effective_offset_ms=%.3f ovr_ms=%.3f\n"
        "  ms_per_frame=%.6f nfrm_in=%d nfrm_out=%d offset_frames=%d effective_offset_frames=%d delta_frames=%d ovr_frames=%d eff_ovr=%d\n"
        "  first_note=%d in_start=%d in_count=%d write_frames=%d keep_out=%d new_nfrm=%d\n",
        data->input, data->output, data->length, data->offset, data->delta,
        effective_offset_ms, data->ovr, ms_per_frame, nfrm_in, nfrm_out,
        offset_frames, effective_offset_frames, delta_frames, ovr_frames, eff_ovr,
        first_note, in_start, in_count, write_frames, keep_out, new_nfrm);

    llsm_container* conf = out_chunk
        ? llsm_copy_container(out_chunk->conf)
        : llsm_copy_container(input->conf);
    llsm_container_attach(conf, LLSM_CONF_NFRM,
        llsm_create_int(new_nfrm), llsm_delete_int, llsm_copy_int);
    // Keep initialized frames for any uncovered gaps in the output timeline.
    llsm_chunk* combined = llsm_create_chunk(conf, 1);
    llsm_delete_container(conf);
    if (!combined) {
        fprintf(stderr, "wavtool: failed to allocate combined chunk\n");
        free(env);
        llsm_delete_chunk(input);
        if (out_chunk) llsm_delete_chunk(out_chunk);
        free(output_llsm);
        free(session_path);
        return 1;
    }

    // Preserve the accumulated phrase. Traditional wavtool mode truncates its
    // overlap tail; phrase mode places requests on an absolute timeline.
    int copy_out = data->use_position ? nfrm_out : keep_out;
    for (int i = 0; i < copy_out; i++) {
        llsm_delete_container(combined->frames[i]);
        combined->frames[i] = llsm_copy_container(out_chunk->frames[i]);
    }


    // Overlap blend + input append
    for (int i = 0; i < in_count; i++) {
        int dst = data->use_position ? absolute_start + i : keep_out + i;
        int out_src = dst;

        int in_overlap = (eff_ovr > 0 && i < eff_ovr && out_src < nfrm_out &&
            combined->frames[dst]);
        if (in_overlap) {
            // alpha ramps 0→1 across the effective overlap window (snapped to output length)
            float alpha = (float)(i + 1) / (float)eff_ovr;
            // The previous note's frames already have its p3/p4 fade-out baked in;
            // avoid applying the new note's attack envelope a second time during
            // the overlap, or the transition dips from double attenuation.
            llsm_container* prev_frame = llsm_copy_container(out_chunk->frames[out_src]);
            llsm_container* env_in = llsm_copy_container(input->frames[in_start + i]);

            apply_gain_db_to_frame(env_in, alias_gain_db);
            ms2_interp_llsm_frame(prev_frame, env_in, alpha);
            float prev_w = 1.0f - alpha;
            float next_w = alpha;
            float norm = sqrtf(prev_w * prev_w + next_w * next_w);
            if (norm > 1e-6f) {
                float comp = 1.0f / norm;
                apply_volume_to_frame(prev_frame, comp * 100.0f);
            }
            llsm_delete_container(combined->frames[dst]);
            combined->frames[dst] = prev_frame;
            llsm_delete_container(env_in);
        } else {
            if (combined->frames[dst]) llsm_delete_container(combined->frames[dst]);
            combined->frames[dst] = llsm_copy_container(input->frames[in_start + i]);
            apply_gain_db_to_frame(combined->frames[dst], alias_gain_db);
            apply_volume_to_frame(combined->frames[dst], env[i]);
        }
    }
    free(env);

    // 6. Save the raw accumulated L1 cache BEFORE phase propagation.
    //    This is what the next wavtool call will load — it must not have
    //    phases baked in or they'll be double-propagated next time.
    llsm_aoptions* opt_a = llsm_create_aoptions();
    FP_TYPE* thop = llsm_container_get(combined->conf, LLSM_CONF_THOP);
    if (thop) opt_a->thop = *thop;
    if (save_llsm(combined, output_llsm, opt_a, &fs, &nbit, 0.0f) != 0)
        fprintf(stderr, "wavtool: failed to save output L1 cache\n");
    llsm_delete_aoptions(opt_a);
    session_append_input(session_path, data->input);
    free(output_llsm);

    if (data->defer_synthesis) {
        llsm_delete_chunk(combined);
        llsm_delete_chunk(input);
        if (out_chunk) llsm_delete_chunk(out_chunk);
        free(session_path);
        free(debug_log);
        return 0;
    }

    // 7. Synthesize the full accumulated phrase from Layer 1
    llsm_soptions* opt_s = llsm_create_soptions((FP_TYPE)fs);
    opt_s->use_l1 = 1;
    llsm_chunk_phasepropagate(combined, 1);
    llsm_output* synth_out = llsm_synthesize(opt_s, combined);
    llsm_delete_soptions(opt_s);

    if (!synth_out || !synth_out->y || ms2_pulse_effects(combined, synth_out)) {
        if (synth_out) llsm_delete_output(synth_out);
        fprintf(stderr, "wavtool: synthesis failed\n");
        append_debug_log(debug_log, 0, "  synth_failed=1\n");
        llsm_delete_chunk(combined);
        llsm_delete_chunk(input);
        if (out_chunk) llsm_delete_chunk(out_chunk);
        free(session_path);
        free(debug_log);
        return 1;
    }

    append_debug_log(debug_log, 0,
        "  synth_samples=%d synth_ms=%.3f expected_phrase_ms=%.3f extra_ms=%.3f\n\n",
        synth_out->ny,
        (double)synth_out->ny * 1000.0 / (double)fs,
        (double)new_nfrm * (double)ms_per_frame,
        (double)synth_out->ny * 1000.0 / (double)fs - (double)new_nfrm * (double)ms_per_frame);

    wavwrite(synth_out->y, synth_out->ny, fs, nbit, data->output);
    llsm_delete_output(synth_out);

    llsm_delete_chunk(combined);
    llsm_delete_chunk(input);
    if (out_chunk) llsm_delete_chunk(out_chunk);
    free(session_path);
    free(debug_log);
    return 0;
}

struct wavtool_phrase {
    llsm_chunk* chunk;
    // Keep fades separate from model parameters until stitching is complete.
    float* envelope;
    int fs;
    int nbit;
};

wavtool_phrase* wavtool_phrase_create(void) {
    return (wavtool_phrase*)calloc(1, sizeof(wavtool_phrase));
}

static llsm_container* silent_frame(const llsm_container* source) {
    llsm_container* frame = llsm_copy_container((llsm_container*)source);
    FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
    if (f0) f0[0] = 0;
    FP_TYPE* vt = llsm_container_get(frame, LLSM_FRAME_VTMAGN);
    if (vt) {
        int count = llsm_fparray_length(vt);
        for (int i = 0; i < count; ++i) vt[i] = -80;
    }
    return frame;
}

int wavtool_phrase_add(wavtool_phrase* phrase, const llsm_chunk* input,
    int fs, int nbit, float gain_db, const wavtool_data* data) {
    double profile_started = ms2_profile_time();
    if (!phrase || !input || !data) return 1;
    int* input_count_ptr = llsm_container_get(input->conf, LLSM_CONF_NFRM);
    if (!input_count_ptr || *input_count_ptr <= 0 || fs <= 0) return 1;
    int input_count = *input_count_ptr;
    float frame_ms = 128000.0f / fs;
    int source_start = max(0, (int)roundf(data->offset / frame_ms));
    if (source_start >= input_count) return 1;
    int take = min(input_count - source_start,
        max(1, (int)roundf(data->length / frame_ms)));
    int destination_start = max(0, (int)roundf(data->position / frame_ms));
    int old_count = 0;
    if (phrase->chunk) {
        int* old_count_ptr = llsm_container_get(phrase->chunk->conf, LLSM_CONF_NFRM);
        old_count = old_count_ptr ? *old_count_ptr : 0;
    }
    int new_count = max(old_count, destination_start + take);
    char profile_detail[128];
    snprintf(profile_detail, sizeof(profile_detail), "old_frames=%d new_frames=%d phone_frames=%d",
        old_count, new_count, take);
    llsm_container* conf = phrase->chunk
        ? llsm_copy_container(phrase->chunk->conf)
        : llsm_copy_container((llsm_container*)input->conf);
    llsm_container_attach(conf, LLSM_CONF_NFRM,
        llsm_create_int(new_count), llsm_delete_int, llsm_copy_int);
    // Every slot is filled below; do not allocate frames only to replace them.
    llsm_chunk* combined = llsm_create_chunk(conf, 0);
    llsm_delete_container(conf);
    if (!combined) return 1;
    for (int i = 0; i < old_count; ++i)
        combined->frames[i] = llsm_copy_container(phrase->chunk->frames[i]);
    for (int i = old_count; i < new_count; ++i)
        combined->frames[i] = silent_frame(input->frames[source_start]);
    ms2_profile("stitch.copy_accumulated_phrase", profile_started, profile_detail);

    float* envelope = build_envelope(data, take, frame_ms);
    if (!envelope) {
        llsm_delete_chunk(combined);
        return 1;
    }
    float* levels = (float*)calloc((size_t)new_count, sizeof(float));
    if (!levels) {
        free(envelope);
        llsm_delete_chunk(combined);
        return 1;
    }
    if (old_count) memcpy(levels, phrase->envelope, sizeof(float) * (size_t)old_count);
    int overlap = max(0, (int)roundf(data->ovr / frame_ms));
    for (int i = 0; i < take; ++i) {
        int destination = destination_start + i;
        llsm_container* next = llsm_copy_container(input->frames[source_start + i]);
        apply_gain_db_to_frame(next, gain_db);
        if (destination < old_count && i < overlap) {
            float previous = levels[destination];
            float incoming = max(0.0f, envelope[i]);
            float total = previous + incoming;
            // The envelope already IS the crossfade. Interpolate unfaded models
            // using its relative weights and apply the summed level only once.
            if (previous <= 0) {
                llsm_delete_container(combined->frames[destination]);
                combined->frames[destination] = next;
            } else {
                if (incoming > 0)
                    ms2_interp_llsm_frame(combined->frames[destination], next, incoming / total);
                llsm_delete_container(next);
            }
            levels[destination] = total;
        } else {
            llsm_delete_container(combined->frames[destination]);
            combined->frames[destination] = next;
            levels[destination] = max(0.0f, envelope[i]);
        }
    }
    free(envelope);
    if (phrase->chunk) llsm_delete_chunk(phrase->chunk);
    free(phrase->envelope);
    phrase->envelope = levels;
    phrase->chunk = combined;
    phrase->fs = fs;
    phrase->nbit = nbit;
    return 0;
}

int wavtool_phrase_synthesize(wavtool_phrase* phrase, float** samples,
    int* sample_count, int* sample_rate) {
    if (!phrase || !phrase->chunk || !samples || !sample_count || !sample_rate) return 1;
    double started = ms2_profile_time();
    llsm_chunk* synthesis = llsm_copy_chunk(phrase->chunk);
    int count = *(int*)llsm_container_get(synthesis->conf, LLSM_CONF_NFRM);
    for (int i = 0; i < count; ++i)
        apply_volume_to_frame(synthesis->frames[i], phrase->envelope[i]);
    llsm_chunk_phasepropagate(synthesis, 1);
    ms2_profile("final.copy_envelopes_phase", started, NULL);
    llsm_soptions* options = llsm_create_soptions((FP_TYPE)phrase->fs);
    options->use_l1 = 1;
    started = ms2_profile_time();
    llsm_output* output = llsm_synthesize(options, synthesis);
    if (output && ms2_pulse_effects(synthesis, output)) {
        llsm_delete_output(output);
        output = NULL;
    }
    ms2_profile("final.llsm_synthesize", started, NULL);
    llsm_delete_soptions(options);
    llsm_delete_chunk(synthesis);
    if (!output || !output->y) {
        if (output) llsm_delete_output(output);
        return 1;
    }
    *samples = (float*)malloc(sizeof(float) * (size_t)output->ny);
    if (!*samples) {
        llsm_delete_output(output);
        return 1;
    }
    memcpy(*samples, output->y, sizeof(float) * (size_t)output->ny);
    *sample_count = output->ny;
    *sample_rate = phrase->fs;
    llsm_delete_output(output);
    return 0;
}

void wavtool_phrase_delete(wavtool_phrase* phrase) {
    if (!phrase) return;
    if (phrase->chunk) llsm_delete_chunk(phrase->chunk);
    free(phrase->envelope);
    free(phrase);
}
