#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include "utils.h"
#include "wavtool.h"
#include "moresampler2.h"
#include <ciglet/ciglet.h>
#include <libllsm/llsm.h>

#ifndef MORESAMPLER2_VERSION
#define MORESAMPLER2_VERSION "dev"
#endif

const char* version = MORESAMPLER2_VERSION;

int base64decoderForUtau(char x, char y)
{
	int ans1, ans2, ans;

	if(x=='+') ans1 = 62;
	if(x=='/') ans1 = 63;
	if(x>='0' && x <= '9') ans1 = x+4;
	if(x>='A' && x <= 'Z') ans1 = x-65;
	if(x>='a' && x <= 'z') ans1 = x-71;

	if(y=='+') ans2 = 62;
	if(y=='/') ans2 = 63;
	if(y>='0' && y <= '9') ans2 = y+4;
	if(y>='A' && y <= 'Z') ans2 = y-65;
	if(y>='a' && y <= 'z') ans2 = y-71;

	ans = (ans1<<6) | ans2;
	if(ans >= 2048) ans -= 4096;
	return ans;
}

int getF0Contour(const char *input, double *output)
{
	int i, j, count, length;
	i = 0;
	count = 0;
	double tmp;
    

	tmp = 0.0;
	while(input[i] != '\0')
	{
		if(input[i] == '#')
		{ 
			length = 0;
			for(j = i+1;input[j]!='#';j++)
			{
				length = length*10 + input[j]-'0';
			}
			i = j+1;
			for(j = 0;j < length;j++)
			{
				output[count++] = tmp;
			}
		}
		else
		{
			tmp = base64decoderForUtau(input[i], input[i+1]);
			output[count++] = tmp;
			i+=2;
		}
	}

	return count;
}

//飴屋／菖蒲氏のworld4utau.cppから移植
double getFreqAvg(double f0[], int tLen)
{
	int i, j;
	double value = 0, r;
	double p[6], q;
	double freq_avg = 0;
	double base_value = 0;
	for (i = 0; i < tLen; i++)
	{
		value = f0[i];
		if (value < 1000.0 && value > 55.0)
		{
			r = 1.0;
			//連続して近い値の場合のウエイトを重くする
			for (j = 0; j <= 5; j++)
			{
				if (i > j) {
					q = f0[i - j - 1] - value;
					p[j] = value / (value + q * q);
				} else {
					p[j] = 1/(1 + value);
				}
				r *= p[j];
			}
			freq_avg += value * r;
			base_value += r;
		}
	}
	if (base_value > 0) freq_avg /= base_value;
	return freq_avg;
}

static int parse_note_to_midi(const char *note_str) {
    // Semitone offsets from C
    int base_note = -1;
    switch (toupper(note_str[0])) {
        case 'C': base_note = 0; break;
        case 'D': base_note = 2; break;
        case 'E': base_note = 4; break;
        case 'F': base_note = 5; break;
        case 'G': base_note = 7; break;
        case 'A': base_note = 9; break;
        case 'B': base_note = 11; break;
        default: return -1;  // invalid note
    }

    int offset = 1;
    if (note_str[offset] == '#') {
        base_note += 1;
        offset++;
    } else if (note_str[offset] == 'b') {
        base_note -= 1;
        offset++;
    }

    int octave = atoi(note_str + offset);
    int midi_note = (octave + 1) * 12 + base_note;
    return midi_note;
}

float note_to_frequency(const char *note_str) {
    int midi = parse_note_to_midi(note_str);
    if (midi < 0) return -1.0f;  // invalid input
    return (float)(440.0 * pow(2.0, (midi - 69) / 12.0));
}

void convert_cents_to_hz_offset(const double* cents, int cents_len,
                                int nfrm, int nhop, int fs, float tempo,
                                float* out_ratio_offset) {

    const float frame_duration_sec = (float)nhop / (float)fs;

    // PIT grid interval from world4utau: pStep samples per PIT point
    const float pit_interval_sec = (60.0f / 96.0f) / tempo;  // seconds per PIT element

    for (int i = 0; i < nfrm; ++i) {
        float time_sec = i * frame_duration_sec;

        float idx = time_sec / pit_interval_sec;

        int i0 = (int)idx;
        if (i0 < 0) i0 = 0;
        if (i0 >= cents_len) i0 = cents_len - 1;

        int i1 = i0 + 1;
        if (i1 >= cents_len) i1 = cents_len - 1;

        float frac = idx - (float)i0;

        float cents_interp = (float)(cents[i0] * (1.0f - frac) + cents[i1] * frac);

        float ratio = powf(2.0f, cents_interp / 1200.0f);
        out_ratio_offset[i] = ratio - 1.0f; // ratio offset, not Hz
    }
}

FP_TYPE* convert_env_to_vol_arr(int* p, int* v, int nfrm) {
    return NULL; // TODO
}

// Port of generate_time_mapping in the original 00413cd0 decompilation:
// output knots [0, fixed * velocity_ratio, duration] map to source knots
// [offset, offset + fixed, end]. Keep the velocity knot even when it lies
// beyond the output: a short note then contains only part of the consonant.
// Coordinates are in frames, but the knots retain fractional frame timing.
static llsm_chunk* map_resampler_frames(llsm_chunk* source, int output_frames,
    double offset, double end, double fixed, double velocity_ratio,
    double duration, int loop) {
  int source_frames = *(int*)llsm_container_get(source->conf, LLSM_CONF_NFRM);
  if (source_frames < 1 || output_frames < 1 || !isfinite(velocity_ratio) ||
      velocity_ratio <= 0 || !isfinite(duration) || duration <= 0) return NULL;
  offset = fmax(0, fmin(offset, source_frames - 1));
  end = fmax(offset, fmin(end, source_frames));
  fixed = fmax(0, fmin(fixed, end - offset));
  double boundary = fixed * velocity_ratio;
  double vowel = end - offset - fixed;

  llsm_container* conf = llsm_copy_container(source->conf);
  llsm_container_attach(conf, LLSM_CONF_NFRM, llsm_create_int(output_frames),
    llsm_delete_int, llsm_copy_int);
  llsm_chunk* out = llsm_create_chunk(conf, 0);
  llsm_delete_container(conf);
  if (!out) return NULL;

  for (int i = 0; i < output_frames; ++i) {
    double position;
    if (i < boundary) {
      position = offset + i / velocity_ratio;
    } else if (vowel <= 0 || duration <= boundary) {
      position = end;
    } else if (loop) {
      // Original loop knots alternate between the vowel start and end.
      double phase = fmod(i - boundary, 2 * vowel);
      position = offset + fixed + (phase <= vowel ? phase : 2 * vowel - phase);
    } else {
      position = offset + fixed + (i - boundary) * vowel / (duration - boundary);
    }
    position = fmax(offset, fmin(position, fmin(end, source_frames - 1)));
    int base = (int)floor(position);
    int next = min(base + 1, source_frames - 1);
    FP_TYPE fraction = (FP_TYPE)(position - base);
    out->frames[i] = llsm_copy_container(source->frames[base]);
    if (!out->frames[i]) {
      llsm_delete_chunk(out);
      return NULL;
    }
    if (next > base && fraction > 0)
      ms2_interp_llsm_frame(out->frames[i], source->frames[next], fraction);
  }
  return out;
}

// Tension is a layer-1 source control. Moresampler's Mt path alters the LF
// source model while rebuilding a frame; libllsm2 exposes the equivalent
// source-shape parameter as Rd. A smaller Rd is tenser/pressed and a larger Rd
// is softer/breathier. Keep a modest tract tilt as the audible spectral part
// of Mt, but do not approximate the entire control with VTMAGN alone.
static void apply_tension_frame(llsm_container* frame, FP_TYPE tension) {
    if (!frame || tension == 0) return;
    if (tension > (FP_TYPE)100.0) tension = (FP_TYPE)100.0;
    if (tension < (FP_TYPE)-100.0) tension = (FP_TYPE)-100.0;
    const FP_TYPE t = tension / (FP_TYPE)100.0;
    const FP_TYPE slope_db = (FP_TYPE)12.0 * t;

    FP_TYPE* rd = (FP_TYPE*)llsm_container_get(frame, LLSM_FRAME_RD);
    if (rd && *rd > 0) {
        FP_TYPE shaped = *rd * (FP_TYPE)pow(2.0, (double)-t);
        if (shaped < (FP_TYPE)0.3) shaped = (FP_TYPE)0.3;
        if (shaped > (FP_TYPE)3.0) shaped = (FP_TYPE)3.0;
        *rd = shaped;
    }

    FP_TYPE* vt = (FP_TYPE*)llsm_container_get(frame, LLSM_FRAME_VTMAGN);
    if (!vt) {
        llsm_container_remove(frame, LLSM_FRAME_HM);
        return;
    }
    int nspec = llsm_fparray_length(vt);
    if (nspec < 2) {
        llsm_container_remove(frame, LLSM_FRAME_HM);
        return;
    }

    // Linear tilt centred at midpoint:
    //   bin 0 (DC)   gets -slope_db/2
    //   bin nspec-1  gets +slope_db/2
    for (int j = 0; j < nspec; ++j) {
        FP_TYPE w = (FP_TYPE)j / (FP_TYPE)(nspec - 1); // 0..1
        vt[j] += slope_db * (w - (FP_TYPE)0.5);
    }
    llsm_container_remove(frame, LLSM_FRAME_HM);
}

static FP_TYPE curve_value(const float* curve, int length, int index, int count) {
    if (!curve || length <= 0 || count <= 1) return 0;
    FP_TYPE position = (FP_TYPE)index * (length - 1) / (count - 1);
    int left = (int)floor(position);
    int right = min(length - 1, left + 1);
    FP_TYPE alpha = position - left;
    return curve[left] * (1 - alpha) + curve[right] * alpha;
}

typedef struct {
    int g, t, Y, e, b, A, P, u;
    int Mt, Mb, Md, Mm, Ms, ME, Mo, Mr, MC, MG, MD, Mp;
} Flags;

void apply_tension(llsm_chunk* chunk, FP_TYPE tension) {
    if (!chunk || tension == 0) return;
    int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
    if (!nfrm_p) return;
    for (int i = 0; i < *nfrm_p; ++i)
        apply_tension_frame(chunk->frames[i], tension);
}

void apply_tension_curve(llsm_chunk* chunk, const float* curve, int length) {
    if (!chunk || !curve || length <= 0) return;
    int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
    if (!nfrm_p) return;
    for (int i = 0; i < *nfrm_p; ++i)
        apply_tension_frame(chunk->frames[i], curve_value(curve, length, i, *nfrm_p));
}


static void apply_gender_frame(llsm_container* frame, FP_TYPE gender) {
  if (!frame || gender == 0) return;
  if (gender > (FP_TYPE)100.0) gender = (FP_TYPE)100.0;
  if (gender < (FP_TYPE)-100.0) gender = (FP_TYPE)-100.0;
  FP_TYPE ratio = (FP_TYPE)pow(2.0, (double)(gender / (FP_TYPE)100.0));


    FP_TYPE* vt = (FP_TYPE*)llsm_container_get(frame, LLSM_FRAME_VTMAGN);
    if (!vt) return;
    int nspec = llsm_fparray_length(vt);
    if (nspec <= 1) return;

    FP_TYPE* src = (FP_TYPE*)malloc(sizeof(FP_TYPE) * (size_t)nspec);
    if (!src) return;
    memcpy(src, vt, sizeof(FP_TYPE) * (size_t)nspec);

    for (int j = 0; j < nspec; ++j) {
      FP_TYPE u = (FP_TYPE)j / (FP_TYPE)(nspec - 1);
      FP_TYPE su = u / ratio;
      if (su < (FP_TYPE)0) su = (FP_TYPE)0;
      if (su > (FP_TYPE)1) su = (FP_TYPE)1;
      FP_TYPE sf = su * (FP_TYPE)(nspec - 1);
      int s0 = (int)floor((double)sf);
      int s1 = min(nspec - 1, s0 + 1);
      FP_TYPE a = sf - (FP_TYPE)s0;
      vt[j] = src[s0] * ((FP_TYPE)1.0 - a) + src[s1] * a;
    }

    free(src);

    llsm_container_remove(frame, LLSM_FRAME_HM);
}

void apply_gender(llsm_chunk* chunk, FP_TYPE gender) {
  if (!chunk || gender == 0) return;
  int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  if (!nfrm_p) return;
  for (int i = 0; i < *nfrm_p; ++i)
    apply_gender_frame(chunk->frames[i], gender);
}

void apply_gender_curve(llsm_chunk* chunk, const float* curve, int length) {
  if (!chunk || !curve || length <= 0) return;
  int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  if (!nfrm_p) return;
  for (int i = 0; i < *nfrm_p; ++i)
    apply_gender_frame(chunk->frames[i], curve_value(curve, length, i, *nfrm_p));
}

static FP_TYPE effect_value(const float* curve, int length, int frame, int count, FP_TYPE fallback) {
  return fallback + (curve && length > 0 ? curve_value(curve, length, frame, count) : 0);
}

static int curve_is_active(const float* curve, int length, float neutral) {
  if (!curve) return 0;
  for (int i = 0; i < length; ++i)
    if (curve[i] != neutral) return 1;
  return 0;
}

static int remaining_effects_active(const ms2_resampler_request* data, const Flags* flags) {
  if (flags->Mb || flags->Md || flags->Mo || flags->ME || flags->Mr ||
      flags->MG || flags->MC || flags->MD || flags->Mp || flags->b ||
      flags->A || flags->Y != 100) return 1;
  const float* curves[] = { data->breathiness_curve, data->mouth_curve,
    data->formant_curve, data->singer_formant_curve, data->growl_curve,
    data->growl_coarseness_curve, data->distortion_curve,
    data->pitch_randomness_curve, data->unvoiced_emphasis_curve,
    data->amplitude_modulation_curve, data->dryness_curve };
  for (size_t i = 0; i < sizeof(curves) / sizeof(curves[0]); ++i)
    if (curve_is_active(curves[i], data->expression_count, 0)) return 1;
  return curve_is_active(data->voicing_curve, data->expression_count, 100);
}

static FP_TYPE clamp_fp(FP_TYPE value, FP_TYPE low, FP_TYPE high) {
  return value < low ? low : (value > high ? high : value);
}

static void apply_noise_gain(llsm_nmframe* noise, FP_TYPE gain_db) {
  if (!noise || gain_db == 0) return;
  for (int i = 0; i < noise->npsd; ++i) noise->psd[i] += gain_db;
  FP_TYPE linear = (FP_TYPE)pow(10.0, (double)(gain_db / 20));
  for (int channel = 0; channel < noise->nchannel; ++channel) {
    noise->edc[channel] *= linear;
    llsm_hmframe* envelope = noise->eenv[channel];
    if (!envelope) continue;
    for (int harmonic = 0; harmonic < envelope->nhar; ++harmonic)
      envelope->ampl[harmonic] *= linear;
  }
}

static void warp_spectrum(FP_TYPE* spectrum, int count, FP_TYPE amount) {
  if (!spectrum || count < 2 || amount == 0) return;
  FP_TYPE ratio = (FP_TYPE)pow(2.0, (double)(amount / 200.0));
  FP_TYPE* source = malloc(sizeof(FP_TYPE) * (size_t)count);
  if (!source) return;
  memcpy(source, spectrum, sizeof(FP_TYPE) * (size_t)count);
  for (int i = 0; i < count; ++i) {
    FP_TYPE position = (FP_TYPE)i / ratio;
    if (position < 0) position = 0;
    if (position > count - 1) position = count - 1;
    int left = (int)position;
    int right = min(count - 1, left + 1);
    FP_TYPE fraction = position - left;
    spectrum[i] = source[left] * (1 - fraction) + source[right] * fraction;
  }
  free(source);
}

static void apply_remaining_effects(llsm_chunk* chunk, const ms2_resampler_request* data,
    const Flags* flags) {
  if (!chunk || !remaining_effects_active(data, flags)) return;
  int* count_ptr = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  FP_TYPE* nyquist_ptr = llsm_container_get(chunk->conf, LLSM_CONF_FNYQ);
  if (!count_ptr) return;
  int count = *count_ptr;
  FP_TYPE nyquist = nyquist_ptr ? *nyquist_ptr : (FP_TYPE)22050.0;
  for (int frame_index = 0; frame_index < count; ++frame_index) {
    llsm_container* frame = chunk->frames[frame_index];
    if (!frame) continue;
    FP_TYPE* vt = llsm_container_get(frame, LLSM_FRAME_VTMAGN);
    llsm_nmframe* noise = llsm_container_get(frame, LLSM_FRAME_NM);
    FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
    int vt_count = vt ? llsm_fparray_length(vt) : 0;
    FP_TYPE breath = effect_value(data->breathiness_curve, data->expression_count,
      frame_index, count, flags->Mb);
    FP_TYPE voicing = data->voicing_curve && data->expression_count > 0
      ? curve_value(data->voicing_curve, data->expression_count, frame_index, count)
      : (FP_TYPE)100.0;
    voicing = clamp_fp(voicing, 0, 100);
    FP_TYPE mouth = effect_value(data->mouth_curve, data->expression_count,
      frame_index, count, flags->Mo);
    FP_TYPE emphasis = effect_value(data->formant_curve, data->expression_count,
      frame_index, count, flags->ME);
    FP_TYPE singer_formant = effect_value(data->singer_formant_curve, data->expression_count,
      frame_index, count, flags->Mr);
    FP_TYPE growl = effect_value(data->growl_curve, data->expression_count,
      frame_index, count, flags->MG);
    FP_TYPE coarse = effect_value(data->growl_coarseness_curve, data->expression_count,
      frame_index, count, flags->MC);
    FP_TYPE distortion = effect_value(data->distortion_curve, data->expression_count,
      frame_index, count, flags->MD);
    FP_TYPE randomness = effect_value(data->pitch_randomness_curve, data->expression_count,
      frame_index, count, flags->Mp);
    FP_TYPE unvoiced_emphasis = effect_value(data->unvoiced_emphasis_curve,
      data->expression_count, frame_index, count, flags->b);
    FP_TYPE amplitude_modulation = effect_value(data->amplitude_modulation_curve,
      data->expression_count, frame_index, count, flags->A);

    warp_spectrum(vt, vt_count, mouth);
    if (vt && vt_count > 2 && emphasis != 0) {
      FP_TYPE* source = malloc(sizeof(FP_TYPE) * (size_t)vt_count);
      if (source) {
        memcpy(source, vt, sizeof(FP_TYPE) * (size_t)vt_count);
        FP_TYPE strength = emphasis / 100;
        for (int i = 1; i < vt_count - 1; ++i) {
          FP_TYPE smooth = (source[i - 1] + source[i] + source[i + 1]) / 3;
          vt[i] += (source[i] - smooth) * strength * 2;
        }
        free(source);
      }
    }
    if (vt && vt_count > 1 && singer_formant != 0) {
      FP_TYPE center = (FP_TYPE)3000.0 / nyquist;
      FP_TYPE gain = singer_formant * (FP_TYPE)0.12;
      for (int i = 0; i < vt_count; ++i) {
        FP_TYPE x = (FP_TYPE)i / (vt_count - 1);
        FP_TYPE distance = (x - center) / (FP_TYPE)0.08;
        vt[i] += gain * (FP_TYPE)exp((double)(-0.5 * distance * distance));
      }
    }
    FP_TYPE dryness = effect_value(data->dryness_curve, data->expression_count,
      frame_index, count, flags->Md);
    FP_TYPE noise_gain = breath * (FP_TYPE)0.12 - dryness * (FP_TYPE)0.12;
    noise_gain += (100 - voicing) * (FP_TYPE)0.18;
    noise_gain += (100 - flags->Y) * (FP_TYPE)0.10;
    if (f0 && *f0 == 0) noise_gain += unvoiced_emphasis * (FP_TYPE)0.12;
    apply_noise_gain(noise, noise_gain);
    FP_TYPE amplitude_gain = 0;
    if (f0 && *f0 > 0 && data->tone > 0 && amplitude_modulation != 0) {
      FP_TYPE pitch_deviation_semitones = 12 * (FP_TYPE)(log(*f0 / data->tone) / log(2.0));
      amplitude_gain = amplitude_modulation * (FP_TYPE)0.12 * pitch_deviation_semitones;
    }
    if (vt) {
      FP_TYPE voiced_gain = -breath * (FP_TYPE)0.06;
      voiced_gain -= max(max(breath, 0), coarse) * (FP_TYPE)0.04;
      FP_TYPE rate = (FP_TYPE)(2.0 + coarse * 0.08);
      FP_TYPE growl_gain = growl * (FP_TYPE)0.08 *
        (FP_TYPE)sin(2.0 * M_PI * rate * frame_index / max(1, count));
      for (int i = 0; i < vt_count; ++i) {
        FP_TYPE normalized = vt_count > 1 ? (FP_TYPE)i / (vt_count - 1) : 0;
        vt[i] += voiced_gain + growl_gain + amplitude_gain;
        vt[i] += coarse * (FP_TYPE)0.035 *
          (FP_TYPE)sin(2.0 * M_PI * (3.0 + coarse * 0.06) * normalized);
        vt[i] += distortion * (FP_TYPE)0.04 *
          (FP_TYPE)sin(2.0 * M_PI * (8.0 + normalized * 24.0) * frame_index / max(1, count));
      }
    }
    if (f0 && *f0 > 0 && randomness > 0) {
      unsigned hash = (unsigned)(frame_index + 1) * 1103515245u + 12345u;
      FP_TYPE noise = ((hash >> 8) & 0xffff) / (FP_TYPE)32767.5 - 1;
      *f0 *= (FP_TYPE)pow(2.0, (double)(noise * randomness * 0.15 / 12.0));
    }
    llsm_container_remove(frame, LLSM_FRAME_HM);
  }
}

// Ms stabilizes analysis artifacts by smoothing layer-1 parameters over time.
// Work from a snapshot so the filter is symmetric and does not drift forward.
static void apply_stabilization(llsm_chunk* chunk, const ms2_resampler_request* data,
    const Flags* flags) {
  if (!chunk || (flags->Ms <= 0 &&
      !curve_is_active(data->stability_curve, data->expression_count, 0))) return;
  int* count_ptr = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  if (!count_ptr || *count_ptr < 2) return;
  int count = *count_ptr;
  llsm_chunk* source = llsm_copy_chunk(chunk);
  if (!source) return;
  for (int frame_index = 0; frame_index < count; ++frame_index) {
    FP_TYPE strength = effect_value(data->stability_curve, data->expression_count,
      frame_index, count, flags->Ms);
    strength = clamp_fp(strength, 0, 10);
    if (strength <= 0) continue;
    int radius = 1 + (int)(strength / 3);
    FP_TYPE alpha = strength / 10;
    llsm_container* dst = chunk->frames[frame_index];
    FP_TYPE* dst_rd = llsm_container_get(dst, LLSM_FRAME_RD);
    FP_TYPE* dst_vt = llsm_container_get(dst, LLSM_FRAME_VTMAGN);
    llsm_nmframe* dst_noise = llsm_container_get(dst, LLSM_FRAME_NM);
    int vt_count = dst_vt ? llsm_fparray_length(dst_vt) : 0;
    FP_TYPE rd_sum = 0;
    int rd_samples = 0;
    FP_TYPE* vt_sum = vt_count > 0 ? calloc((size_t)vt_count, sizeof(FP_TYPE)) : NULL;
    FP_TYPE* noise_sum = dst_noise && dst_noise->npsd > 0
      ? calloc((size_t)dst_noise->npsd, sizeof(FP_TYPE)) : NULL;
    int spectral_samples = 0;
    for (int nearby = max(0, frame_index - radius);
         nearby <= min(count - 1, frame_index + radius); ++nearby) {
      llsm_container* src = source->frames[nearby];
      FP_TYPE* src_rd = llsm_container_get(src, LLSM_FRAME_RD);
      FP_TYPE* src_vt = llsm_container_get(src, LLSM_FRAME_VTMAGN);
      llsm_nmframe* src_noise = llsm_container_get(src, LLSM_FRAME_NM);
      if (dst_rd && src_rd && *src_rd > 0) { rd_sum += *src_rd; ++rd_samples; }
      if (vt_sum && src_vt && llsm_fparray_length(src_vt) == vt_count &&
          noise_sum && src_noise && src_noise->npsd == dst_noise->npsd) {
        for (int i = 0; i < vt_count; ++i) vt_sum[i] += src_vt[i];
        for (int i = 0; i < dst_noise->npsd; ++i) noise_sum[i] += src_noise->psd[i];
        ++spectral_samples;
      }
    }
    if (dst_rd && rd_samples > 0)
      *dst_rd = *dst_rd * (1 - alpha) + rd_sum / rd_samples * alpha;
    if (spectral_samples > 0) {
      for (int i = 0; i < vt_count; ++i)
        dst_vt[i] = dst_vt[i] * (1 - alpha) + vt_sum[i] / spectral_samples * alpha;
      for (int i = 0; i < dst_noise->npsd; ++i)
        dst_noise->psd[i] = dst_noise->psd[i] * (1 - alpha) +
          noise_sum[i] / spectral_samples * alpha;
    }
    free(vt_sum);
    free(noise_sum);
    llsm_container_remove(dst, LLSM_FRAME_HM);
  }
  llsm_delete_chunk(source);
}

FP_TYPE* get_pitch_from_area(llsm_chunk* chunk, int start, int end) {
  int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  if (!nfrm_p) return NULL;

  FP_TYPE* pitch = malloc(sizeof(FP_TYPE) * (end - start));
  if (!pitch) return NULL;

  for (int i = start; i < end; ++i) {
    FP_TYPE f0 = *((FP_TYPE*)llsm_container_get(chunk->frames[i], LLSM_FRAME_F0));
    pitch[i - start] = f0;
  }

  return pitch;
}

FP_TYPE* interp_pitch_to_pitd(FP_TYPE* pitch, int old_nfrm, int new_nfrm) {
  FP_TYPE* new_pitch = malloc(sizeof(FP_TYPE) * new_nfrm);
  if (!new_pitch) return NULL;

  for (int i = 0; i < new_nfrm; ++i) {
    float idx = (float)i * (float)(old_nfrm - 1) / (float)(new_nfrm - 1);
    int idx0 = (int)idx;
    int idx1 = min(idx0 + 1, old_nfrm - 1);
    float frac = idx - (float)idx0;
    new_pitch[i] = pitch[idx0] * (1.0f - frac) + pitch[idx1] * frac;
  }
  
  // Get median ignoring 0 values
  FP_TYPE* nonzero = malloc(sizeof(FP_TYPE) * new_nfrm);
  int nz_count = 0;
  for (int i = 0; i < new_nfrm; ++i) {
    if (new_pitch[i] > 0.0f) {
      nonzero[nz_count++] = new_pitch[i];
    }
  }

  FP_TYPE average = (nz_count > 0) ? nonzero[nz_count / 2] : 0.0f;

  free(nonzero);

  for (int i = 0; i < new_nfrm; ++i) {
    new_pitch[i] = new_pitch[i] - average;
  }

  return new_pitch;
}
//modulation of original sample's pitch contour
FP_TYPE* apply_modulation(FP_TYPE* pitch, int mod, int nfrm) {
  FP_TYPE* modulated_pitch = malloc(sizeof(FP_TYPE) * nfrm);
  if (!modulated_pitch) return NULL;

  for (int i = 0; i < nfrm; ++i) {
    modulated_pitch[i] = pitch[i] * (mod / 100.0f);
  }
  return modulated_pitch;
}

int clamp(int val, int min, int max) {
    return val < min ? min : (val > max ? max : val);
}

void parse_flag_string(const char* str, Flags* flags_out) {
    memset(flags_out, 0, sizeof(*flags_out));
    // Optional normalization is opt-in, just like the other effects.
    flags_out->P = 0;
    flags_out->Y = 100;
    flags_out->Mm = 100;

    const char* p = str;
    while (*p != '\0') {
        char* end;
        if (strncmp(p, "Me", 2) == 0) {
            p += 2;
            flags_out->e = 1;
        } else if (strncmp(p, "Mt", 2) == 0 || strncmp(p, "Mb", 2) == 0 ||
            strncmp(p, "Md", 2) == 0 || strncmp(p, "Mm", 2) == 0 ||
            strncmp(p, "Ms", 2) == 0 || strncmp(p, "ME", 2) == 0 ||
            strncmp(p, "Mo", 2) == 0 || strncmp(p, "Mr", 2) == 0 ||
            strncmp(p, "MC", 2) == 0 || strncmp(p, "MG", 2) == 0 ||
            strncmp(p, "MD", 2) == 0 || strncmp(p, "Mp", 2) == 0) {
            char a = p[1]; p += 2;
            int value = (int)strtol(p, &end, 10); p = end;
            if (a == 't') flags_out->Mt = clamp(value, -100, 100);
            else if (a == 'b') flags_out->Mb = clamp(value, -100, 100);
            else if (a == 'd') flags_out->Md = clamp(value, -100, 100);
            else if (a == 'm') flags_out->Mm = clamp(value, 0, 100);
            else if (a == 's') flags_out->Ms = clamp(value, 0, 10);
            else if (a == 'E') flags_out->ME = clamp(value, -100, 100);
            else if (a == 'o') flags_out->Mo = clamp(value, -100, 100);
            else if (a == 'r') flags_out->Mr = clamp(value, -100, 100);
            else if (a == 'C') flags_out->MC = clamp(value, 0, 100);
            else if (a == 'G') flags_out->MG = clamp(value, 0, 100);
            else if (a == 'D') flags_out->MD = clamp(value, 0, 100);
            else if (a == 'p') flags_out->Mp = clamp(value, 0, 100);
        } else if (*p == 't') {
            p++;
            flags_out->t = strtol(p, &end, 10);
            p = end;
            flags_out->t = clamp(flags_out->t, -1200, 1200);
        } else if (*p == 'P') {
            p++;
            flags_out->P = strtol(p, &end, 10);
            p = end;
            flags_out->P = clamp(flags_out->P, 0, 100);
        } else if (*p == 'g') {
            p++;
            flags_out->g = strtol(p, &end, 10);
            p = end;
            flags_out->g = clamp(flags_out->g, -100, 100);
        } else if (*p == 'e') {
            p++;
            flags_out->e = 0;
        } else if (*p == 'b' || *p == 'A' || *p == 'Y') {
            char name = *p++; int value = (int)strtol(p, &end, 10); p = end;
            if (name == 'b') flags_out->b = clamp(value, -20, 100);
            else if (name == 'A') flags_out->A = clamp(value, -100, 100);
            else flags_out->Y = clamp(value, 0, 100);
        } else if (*p == 'u') {
            p++; flags_out->u = 1;
        } else {
            p++;
        }
    }
}

int parse_tempo(const char *tempo_str) {
    if (tempo_str[0] == '!') {
        tempo_str++;  // skip the '!'
    }
    return atoi(tempo_str);
}

static int ms2_process(const ms2_resampler_request* data,
    llsm_chunk** phrase_chunk, float* phrase_gain_db,
    int* phrase_fs, int* phrase_nbit) {
  double profile_started = ms2_profile_time();
  ms2_profile("phone.begin", profile_started, data->input);
  // Allocate and load pitch curve
  int direct_pitch_count = data->pitch_values && data->pitch_value_count > 0
    ? data->pitch_value_count : 0;
  int f0_capacity = direct_pitch_count > 0 ? direct_pitch_count : 3000;
  double* f0_curve = malloc(sizeof(double) * (size_t)f0_capacity);
  if (!f0_curve) return 1;
  int pit_len = direct_pitch_count;
  if (direct_pitch_count > 0) {
    for (int i = 0; i < direct_pitch_count; ++i) f0_curve[i] = data->pitch_values[i];
  } else if (data->pitch_curve && data->pitch_curve[0]) {
    pit_len = getF0Contour(data->pitch_curve, f0_curve);
  } else {
    f0_curve[0] = 0;
    pit_len = 1;
  }
  if (!pit_len) { free(f0_curve); return 1; }

  float velocity = (float)exp2(1 - data->velocity / 100.0f);
  Flags flags;
  parse_flag_string(data->flags, &flags);

  llsm_aoptions* opt_a = llsm_create_aoptions();
  llsm_soptions* opt_s = llsm_create_soptions(44100.0f);

  int nhop = 128;
  int fs = 0, nbit = 0, nx = 0;
  int nfrm = 0;

  double load_started = ms2_profile_time();
  llsm_chunk* chunk;
  if (data->analysis_file && data->analysis_file[0]) {
    chunk = read_llsm_options(data->analysis_file, &nfrm, &fs, &nbit, NULL, opt_a);
    ms2_profile(chunk ? "analysis.cache_hit" : "analysis.cache_invalid", load_started, data->analysis_file);
    if (chunk && !ms2_chunk_is_layer1(chunk)) { llsm_delete_chunk(chunk); chunk = NULL; }
  } else {
    chunk = get_chunk_from_file(data->input, &nfrm, &fs, &nbit,
      &nx, &nhop, 0, 0, opt_a);
  }
  if (!chunk) {
    free(f0_curve);
    return 1;
  }

  printf("Phase sync/stretching\n");
  // Calculate start and end frames based on offset and cutoff (in ms)
  int start_frame = (int)round((data->offset / 1000.0) * fs / nhop);
  int end_frame;
  if (data->cutoff < 0) {
      // Negative cutoff: measured from offset
      end_frame = (int)round(((data->offset + fabs(data->cutoff)) / 1000.0) * fs / nhop);
  } else {
      // Positive cutoff: measured from end of file
      end_frame = nfrm - (int)round((data->cutoff / 1000.0) * fs / nhop);
  }
  if (start_frame < 0) start_frame = 0;
  if (start_frame >= nfrm) start_frame = nfrm - 1;
  if (end_frame > nfrm) end_frame = nfrm;
  if (end_frame <= start_frame) end_frame = start_frame + 1;

  
  int sample_frames = end_frame - start_frame;
  double frames_per_ms = fs / (1000.0 * nhop);
  double output_duration = data->length * frames_per_ms;
  int total_frames = (int)round(output_duration);
  if (total_frames < 1) total_frames = 1;
  if (output_duration <= 0) output_duration = 1;

  float* f0_array = malloc(sizeof(float) * total_frames);
  if (!f0_array) {
    // Handle out-of-memory error
    free(f0_curve);
    if (chunk) llsm_delete_chunk(chunk);
    if (opt_a) llsm_delete_aoptions(opt_a);
    if (opt_s) llsm_delete_soptions(opt_s);
    return 1;
  }

  FP_TYPE* mod_pitch = data->modulation == 0
    ? calloc((size_t)total_frames, sizeof(FP_TYPE))
    : apply_modulation(interp_pitch_to_pitd(get_pitch_from_area(chunk, start_frame, end_frame), sample_frames, total_frames), data->modulation, total_frames);
  
  convert_cents_to_hz_offset(f0_curve, pit_len, total_frames, nhop, fs, data->tempo, f0_array);
  printf("pit_len: %d, total_frames: %d\n", pit_len, total_frames);
  // Convert and synchronize the complete source before mapping. Truncating it
  // to output length first loses frames needed by fast consonants/short notes.
  double stage_started = ms2_profile_time();
  // Source caches are already L1; do not repeat glottal/spectral analysis.
  llsm_chunk_phasepropagate(chunk, -1);
  ms2_profile("phone.source_layer1_phase", stage_started, data->input);
  stage_started = ms2_profile_time();
  double source_end = data->cutoff < 0
    ? (data->offset - data->cutoff) * frames_per_ms
    : nfrm - data->cutoff * frames_per_ms;
  llsm_chunk* chunk_new = map_resampler_frames(chunk, total_frames,
    data->offset * frames_per_ms, source_end,
    data->consonant * frames_per_ms, velocity, output_duration, flags.e == 1);
  if (!chunk_new) {
    free(f0_array);
    free(mod_pitch);
    free(f0_curve);
    llsm_delete_chunk(chunk);
    llsm_delete_aoptions(opt_a);
    llsm_delete_soptions(opt_s);
    return 1;
  }
  ms2_profile("phone.time_mapping", stage_started, data->input);
  stage_started = ms2_profile_time();
  // Set all f0 to target tone
  for(int i = 0; i < total_frames; i ++) {
    llsm_container_attach(chunk_new->frames[i], LLSM_FRAME_HM, NULL, NULL, NULL);
    FP_TYPE* f0_i = llsm_container_get(chunk_new->frames[i], LLSM_FRAME_F0);
    FP_TYPE old_f0 = f0_i[0];

    if (f0_i[0] == 0.00f) {
      continue;
    }
    else {
      f0_i[0] = data->tone * pow(2.0f, (double)flags.t / 1200.0) * (1.0 + f0_array[i]);
      f0_i[0] += mod_pitch[i];
      if (f0_i[0] < 20.0f && f0_i[0] != 0.0f) f0_i[0] = 20.0f; // minimum F0
    }
    // Compensate for the amplitude gain.
    FP_TYPE* vt_magn = llsm_container_get(chunk_new->frames[i],
      LLSM_FRAME_VTMAGN);
    if(vt_magn != NULL) {
      int nspec = llsm_fparray_length(vt_magn);
      for(int j = 0; j < nspec; j ++)
        vt_magn[j] -= 20.0 * log10(f0_i[0] / old_f0);
    }
    
  }
  free(f0_array);
  free(mod_pitch);
  apply_gender(chunk_new, flags.g);
  if (curve_is_active(data->gender_curve, data->expression_count, 0))
    apply_gender_curve(chunk_new, data->gender_curve, data->expression_count);
  apply_tension(chunk_new, flags.Mt);
  if (curve_is_active(data->tension_curve, data->expression_count, 0))
    apply_tension_curve(chunk_new, data->tension_curve, data->expression_count);
  apply_remaining_effects(chunk_new, data, &flags);
  apply_stabilization(chunk_new, data, &flags);
  ms2_profile("phone.pitch_effects", stage_started, data->input);
  stage_started = ms2_profile_time();
  // Phrase mode hands L1 frames directly to the stitcher. Only explicit peak
  // normalization needs temporary PCM to measure the per-phone peak, and that
  // PCM stays in memory. Standalone resampling still produces its usual files.
  int synthesize_phone = !phrase_chunk || flags.P > 0;
  llsm_chunk* cache_chunk = synthesize_phone ? llsm_copy_chunk(chunk_new) : chunk_new;
  llsm_output* out = NULL;
  float alias_gain_db = 0;
  if (synthesize_phone) {
    llsm_chunk_tolayer0(chunk_new);
    llsm_chunk_phasepropagate(chunk_new, 1);
    printf("Synthesis\n");
    out = llsm_synthesize(opt_s, chunk_new);
    if (!out || !out->y) {
      printf("Failed to synthesize output\n");
      if (cache_chunk) llsm_delete_chunk(cache_chunk);
      if (out) llsm_delete_output(out);
      llsm_delete_chunk(chunk_new);
      llsm_delete_chunk(chunk);
      llsm_delete_aoptions(opt_a);
      llsm_delete_soptions(opt_s);
      free(f0_curve);
      return 1;
    }
    alias_gain_db = normalize_waveform(out->y, out->ny, 0.50f, flags.P);
  }

  ms2_profile(synthesize_phone ? "phone.pcm_synthesis" : "phone.pcm_skipped", stage_started, data->input);
  if (phrase_chunk) {
    float volume = data->volume > 0 ? data->volume / 100.0f : 0.0001f;
    *phrase_chunk = cache_chunk;
    if (phrase_gain_db) *phrase_gain_db = alias_gain_db + 20.0f * log10f(volume);
    if (phrase_fs) *phrase_fs = fs;
    if (phrase_nbit) *phrase_nbit = nbit;
    if (out) llsm_delete_output(out);
    llsm_delete_chunk(chunk);
    if (chunk_new != cache_chunk) llsm_delete_chunk(chunk_new);
    llsm_delete_aoptions(opt_a);
    llsm_delete_soptions(opt_s);
    free(f0_curve);
    ms2_profile("phone.total", profile_started, data->input);
    return 0;
  }

  char* output_llsm_path = build_llsm_path(data->output);
  printf("Saving synthesis result to cache: %s\n", output_llsm_path);
  if (!cache_chunk || save_llsm(cache_chunk, output_llsm_path, opt_a,
      &fs, &nbit, alias_gain_db) != 0) {
    printf("Failed to save chunk.\n");
  }
  free(output_llsm_path);
  if (cache_chunk) llsm_delete_chunk(cache_chunk);

  float scale = data->volume / 100.0f;
  for (int i = 0; i < out->ny; ++i)
    out->y[i] *= scale;
  
  wavwrite(out->y, out->ny, fs, nbit, (char*)data->output);
  

  llsm_delete_output(out);
  llsm_delete_chunk(chunk);
  llsm_delete_chunk(chunk_new);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0_curve);
  return 0;
}

int ms2_resample(const ms2_resampler_request* data) {
  return ms2_process(data, NULL, NULL, NULL, NULL);
}

int ms2_analyze_v1(const char* input, const char* output) {
  return ms2_analyze_config_v1(input, output, NULL);
}

int ms2_analyze_config_v1(const char* input, const char* output,
    const ms2_analysis_settings* settings) {
  return ms2_generate_v1(input, output, settings, 0);
}

int ms2_generate_v1(const char* input, const char* output,
    const ms2_analysis_settings* settings, int overwrite) {
  if (!input || !input[0] || !output || !output[0]) return 1;
  struct stat cache_stat;
  if (!overwrite && stat(output, &cache_stat) == 0) {
    if (ms2_validate_analysis_v1(output) == 0) return 0;
    fprintf(stderr, "WARNING: Unusable analysis; regenerating %s\n", output);
    ms2_profile("analysis.regenerate_invalid", ms2_profile_time(), output);
  }
  if (settings && (!isfinite(settings->f0_min) || !isfinite(settings->f0_max) ||
      !isfinite(settings->voicing_bias) || settings->f0_min < 20 || settings->f0_max > 2000 ||
      settings->f0_min >= settings->f0_max || settings->voicing_bias < .1f ||
      settings->voicing_bias > 10 || settings->transition_range < 1 ||
      settings->transition_range > 48 || (settings->refine_f0 != 0 && settings->refine_f0 != 1))) return 4;
  llsm_aoptions* options = llsm_create_aoptions();
  int frames = 0, fs = 0, bits = 0, samples = 0, hop = 128;
  llsm_chunk* chunk = get_chunk_config(input, &frames, &fs, &bits,
    &samples, &hop, 1, 1, options, settings);
  if (chunk && ms2_prepare_layer1(chunk)) { llsm_delete_chunk(chunk); chunk = NULL; }
  if (!chunk) {
    llsm_delete_aoptions(options);
    return 2;
  }
  int code = ms2_save_analysis(chunk, output, options, &fs, &bits) == 0 ? 0 : 3;
  llsm_delete_chunk(chunk);
  llsm_delete_aoptions(options);
  return code;
}

int ms2_validate_analysis_v1(const char* analysis) {
  double started = ms2_profile_time();
  if (!analysis || !analysis[0]) return 1;
  int frames = 0, fs = 0, bits = 0;
  llsm_chunk* chunk = read_llsm(analysis, &frames, &fs, &bits, NULL);
  ms2_profile(chunk ? "analysis.validate_loaded" : "analysis.validate_failed", started, analysis);
  if (!chunk) return 2;
  int valid = frames > 0 && fs > 0 && ms2_chunk_is_layer1(chunk);
  llsm_delete_chunk(chunk);
  return valid ? 0 : 3;
}

int ms2_render_phrase_v1(const ms2_phrase_request* request, ms2_render_result* result) {
  if (!request || !result || !request->phones || request->phone_count <= 0 ||
      !request->output || !request->output[0]) return 1;
  result->samples = NULL;
  result->sample_count = 0;
  result->sample_rate = 0;
  double profile_started = ms2_profile_time();
  ms2_profile("phrase.begin", profile_started, request->output);
  wavtool_phrase* phrase = wavtool_phrase_create();
  if (!phrase) return 1;

  for (int p = 0; p < request->phone_count; ++p) {
    const ms2_phrase_phone* phone = &request->phones[p];
    llsm_chunk* phone_chunk = NULL;
    float gain_db = 0;
    int fs = 0, nbit = 0;
    if (ms2_process(&phone->synthesis, &phone_chunk, &gain_db, &fs, &nbit) != 0) {
      wavtool_phrase_delete(phrase);
      return 2;
    }
    if (!phone->envelope_x_ms || !phone->envelope_y || phone->envelope_count < 5) {
      llsm_delete_chunk(phone_chunk);
      wavtool_phrase_delete(phrase);
      return 3;
    }
    wavtool_data data = {0};
    data.offset = phone->skip_ms;
    data.length = phone->length_ms;
    data.p1 = 0;
    data.p2 = phone->envelope_x_ms[1] - phone->envelope_x_ms[0];
    data.p3 = phone->envelope_x_ms[4] - phone->envelope_x_ms[3];
    data.p4 = 0;
    data.p5 = phone->envelope_x_ms[2] - phone->envelope_x_ms[0];
    data.v1 = phone->envelope_y[0] * 100;
    data.v2 = phone->envelope_y[1] * 100;
    data.v3 = phone->envelope_y[3] * 100;
    data.v4 = phone->envelope_y[4] * 100;
    data.v5 = phone->envelope_y[2] * 100;
    data.ovr = phone->overlap_ms;
    data.position = phone->position_ms;
    data.use_position = 1;
    double stitch_started = ms2_profile_time();
    if (wavtool_phrase_add(phrase, phone_chunk, fs, nbit, gain_db, &data) != 0) {
      llsm_delete_chunk(phone_chunk);
      wavtool_phrase_delete(phrase);
      return 4;
    }
    ms2_profile("phrase.stitch_phone", stitch_started, request->output);
    llsm_delete_chunk(phone_chunk);
  }

  int rendered_count = 0;
  float* rendered = NULL;
  double synth_started = ms2_profile_time();
  if (wavtool_phrase_synthesize(phrase, &rendered, &rendered_count,
      &result->sample_rate) != 0) {
    wavtool_phrase_delete(phrase);
    return 5;
  }
  wavtool_phrase_delete(phrase);
  ms2_profile("phrase.final_synthesis", synth_started, request->output);
  int target_count = (int)ceil(request->duration_ms * result->sample_rate / 1000.0f);
  if (target_count <= 0) {
    free(rendered);
    return 5;
  }
  result->samples = (float*)calloc((size_t)target_count, sizeof(float));
  if (!result->samples) {
    free(rendered);
    return 5;
  }
  memcpy(result->samples, rendered,
    sizeof(float) * (size_t)min(rendered_count, target_count));
  free(rendered);
  result->sample_count = target_count;
  double write_started = ms2_profile_time();
  wavwrite(result->samples, result->sample_count, result->sample_rate, 16,
    (char*)request->output);
  ms2_profile("phrase.write_wav", write_started, request->output);
  ms2_profile("phrase.total", profile_started, request->output);
  return 0;
}

void ms2_free(void* pointer) {
  free(pointer);
}

#ifndef MORESAMPLER2_NO_MAIN
int main(int argc, char* argv[]) {
    printf("moresampler2 version %s\n", version);
    if (argc == 2) { // user dragged and dropped a folder into the executable
        printf("At the moment, autolabeling is not supported.\n");
        return 0;
    }
    if (argc < 2) {
        printf("Moresampler is meant to be used inside of UTAU or OpenUtau.\n");
        return 1;
    }
    // Wavtool mode: argc 14–16, argv[3] is a numeric offset (ms).
    // Resampler mode: argc 14, argv[3] is a note name (e.g. "A4").
    // When argc == 14, disambiguate by checking if argv[3] starts with a digit.
    if (argc >= 14 && argc <= 16 &&
        (isdigit((unsigned char)argv[3][0]) || argv[3][0] == '-')) {
        wavtool_data data = {0};
        data.output = argv[1];
        data.input  = argv[2];
        data.offset = atof(argv[3]);
        data.length = atof(argv[4]);
        data.p1     = atof(argv[5]);
        data.p2     = atof(argv[6]);
        data.p3     = atof(argv[7]);
        data.v1     = atof(argv[8]);
        data.v2     = atof(argv[9]);
        data.v3     = atof(argv[10]);
        data.v4     = atof(argv[11]);
        data.ovr    = atof(argv[12]);
        data.p4     = atof(argv[13]);
        data.p5     = (argc >= 15) ? atof(argv[14]) : 0.0f;
        data.v5     = (argc >= 16) ? atof(argv[15]) : 0.0f;
        return wavtool(&data);
    }
    if (argc == 14) { // user wants the resampler mode
        ms2_resampler_request data = {0};
        data.input = argv[1];
        data.output = argv[2];
        data.tone = note_to_frequency(argv[3]); // note is passed as a string (A4), so we need to convert it
        data.velocity = atof(argv[4]);
        data.flags = argv[5]; // flags are passed as a string, e.g. "Mt50"
        data.offset = atof(argv[6]);
        data.length = atof(argv[7]);
        data.consonant = atof(argv[8]);
        data.cutoff = atof(argv[9]);
        data.volume = atoi(argv[10]);
        data.modulation = atoi(argv[11]);
        data.tempo = parse_tempo(argv[12]); //since tempo has a special format, we need to parse it
        data.pitch_curve = argv[13]; // pitch curve data as a string
        return ms2_resample(&data);
    }
    printf("Invalid arguments. Expected 14 arguments, got %d.\n", argc);
    return 0;
}
#endif
