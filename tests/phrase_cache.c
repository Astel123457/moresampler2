#include <libllsm/llsm.h>
#ifdef _WIN32
#include <sys/utime.h>
#else
#include <utime.h>
#endif
static int synthesis_calls;
static llsm_output* counted_synthesis(llsm_soptions* options, llsm_chunk* chunk) {
    ++synthesis_calls;
    return llsm_synthesize(options, chunk);
}
#define llsm_synthesize counted_synthesis
#define MORESAMPLER2_NO_MAIN
#include "../moresampler2.c"
#undef llsm_synthesize

static void check(int ok, const char* message) {
    if (!ok) { fprintf(stderr, "%s\n", message); exit(1); }
}

static int exists(const char* path) {
    FILE* file = fopen(path, "rb");
    if (!file) return 0;
    fclose(file);
    return 1;
}

int main(void) {
    const char* profile = "phrase-cache-profile.log";
    remove(profile);
    ms2_set_profile_path_v1(profile);
    const char* input = "phrase-cache-source.wav";
    const char* analysis = "phrase-cache-source.llsm2";
    const char* intermediate = "phrase-cache-intermediate.wav";
    const char* intermediate_model = "phrase-cache-intermediate.llsm2";
    remove(analysis);
    remove(intermediate);
    remove(intermediate_model);
    FP_TYPE wave[22050];
    for (int i = 0; i < 22050; ++i)
        wave[i] = .2f * sin(2 * M_PI * 220 * i / 44100.0);
    wavwrite(wave, 22050, 44100, 16, (char*)input);
    ms2_resampler_request request = {0};
    request.input = input;
    request.output = intermediate;
    request.tone = 220;
    request.velocity = request.volume = 100;
    request.length = 200;
    request.consonant = 40;
    request.cutoff = -250;
    request.tempo = 120;
    request.flags = "";
    for (int pass = 0; pass < 3; ++pass) {
        // Last pass explicitly requests normalization: PCM may exist in memory,
        // but must still never be written as intermediate WAV/model files.
        if (pass == 2) request.flags = "P86";
        llsm_chunk* chunk = NULL;
        float gain = 0;
        int fs = 0, bits = 0;
        check(!ms2_process(&request, &chunk, &gain, &fs, &bits), "phrase phone failed");
        check(chunk && fs == 44100, "phrase phone must return model frames");
        check(synthesis_calls == (pass == 2 ? 1 : 0), "P0 must bypass per-phone synthesis");
        check(exists(analysis), "source analysis must remain on disk");
        check(!exists(intermediate) && !exists(intermediate_model), "phrase mode wrote an intermediate file");
        llsm_delete_chunk(chunk);
    }
    check(ms2_validate_analysis_v1(analysis) == 0, "generated analysis must validate as Layer 1");
    // Exercise repeated cache loading without re-analysis. With allocation
    // accounting enabled, any lost decoded/preinitialized frames fail at exit.
    for (int pass = 0; pass < 32; ++pass) {
        check(ms2_validate_analysis_v1(analysis) == 0, "repeated cache validation failed");
    }
    int count = 0, fs = 0, bits = 0;
    llsm_aoptions* options = llsm_create_aoptions();
    llsm_chunk* legacy = read_llsm_options(analysis, &count, &fs, &bits, NULL, options);
    check(legacy && ms2_chunk_is_layer1(legacy), "saved source must contain Layer 1");
    llsm_chunk_tolayer0(legacy);
    for (int i = 0; i < count; ++i) {
        if (!llsm_container_get(legacy->frames[i], LLSM_FRAME_HM))
            llsm_container_attach(legacy->frames[i], LLSM_FRAME_HM,
                llsm_create_hmframe(0), llsm_delete_hmframe, llsm_copy_hmframe);
        llsm_container_remove(legacy->frames[i], LLSM_FRAME_RD);
        llsm_container_remove(legacy->frames[i], LLSM_FRAME_VTMAGN);
        llsm_container_remove(legacy->frames[i], LLSM_FRAME_VSPHSE);
    }
    check(!save_llsm(legacy, analysis, options, &fs, &bits, 0), "legacy fixture save failed");
    llsm_delete_chunk(legacy);
    llsm_delete_aoptions(options);
    check(ms2_validate_analysis_v1(analysis) != 0, "Singers validator must identify legacy Layer 0");
    check(!ms2_analyze_v1(input, analysis), "legacy analysis must regenerate automatically");
    check(!ms2_validate_analysis_v1(analysis), "automatic generation must upgrade legacy analysis");
    check(!ms2_generate_v1(input, analysis, NULL, 1), "explicit overwrite must regenerate analysis");
    check(!ms2_validate_analysis_v1(analysis), "upgraded analysis must validate");
    check(!ms2_analyze_v1(input, analysis), "repeated generation must reuse Layer 1");
    // Non-default generation must honor native options and never overwrite the
    // default source cache. Renderer requests explicitly select the variant.
    const char* variant = "phrase-cache-configured.llsm2";
    remove(variant);
    ms2_analysis_settings configured = {100, 400, 1.5f, 12, 0};
    check(!ms2_analyze_config_v1(input, variant, &configured), "configured analysis failed");
    options = llsm_create_aoptions();
    legacy = read_llsm_options(variant, &count, &fs, &bits, NULL, options);
    check(legacy && ms2_chunk_is_layer1(legacy) && options->f0_refine == 0,
        "configured analysis must store Layer 1 with refinement disabled");
    llsm_delete_chunk(legacy);
    llsm_delete_aoptions(options);
    configured.f0_min = 500;
    check(ms2_generate_v1(input, variant, &configured, 1) == 4, "invalid F0 interval must fail when generating");
    configured.f0_min = 120;
    configured.refine_f0 = 1;
#ifdef _WIN32
    struct _utimbuf old_time = {1, 1};
    check(_utime(variant, &old_time) == 0, "could not mark analysis older than source");
#else
    struct utimbuf old_time = {1, 1};
    check(utime(variant, &old_time) == 0, "could not mark analysis older than source");
#endif
    // Reuse must not touch the file, even when the source is newer and settings differ.
    FILE* source_again = fopen(input, "ab");
    check(source_again != NULL, "source reopen failed");
    fputc(0, source_again);
    fclose(source_again);
    check(!ms2_analyze_config_v1(input, variant, &configured), "existing analysis must be reused");
    options = llsm_create_aoptions();
    legacy = read_llsm_options(variant, &count, &fs, &bits, NULL, options);
    check(legacy && options->f0_refine == 0, "settings changes must not replace existing analysis");
    llsm_delete_chunk(legacy);
    llsm_delete_aoptions(options);
    check(!ms2_generate_v1(input, variant, &configured, 1), "manual overwrite must apply new settings");
    options = llsm_create_aoptions();
    legacy = read_llsm_options(variant, &count, &fs, &bits, NULL, options);
    check(legacy && options->f0_refine == 1, "manual overwrite ignored changed settings");
    llsm_delete_chunk(legacy);
    llsm_delete_aoptions(options);
    remove(analysis);
    request.analysis_file = variant;
    request.flags = "";
    llsm_chunk* configured_phone = NULL;
    check(!ms2_process(&request, &configured_phone, NULL, NULL, NULL), "renderer must read selected analysis variant");
    check(!exists(analysis), "renderer must not regenerate the default analysis when using a variant");
    llsm_delete_chunk(configured_phone);
    FILE* broken = fopen(variant, "wb");
    check(broken != NULL, "could not create invalid analysis fixture");
    fputs("invalid", broken);
    fclose(broken);
    check(!ms2_analyze_config_v1(input, variant, &configured),
        "invalid existing analysis must regenerate automatically");
    check(!ms2_validate_analysis_v1(variant), "regenerated analysis must validate");
    broken = fopen(variant, "wb");
    check(broken != NULL, "could not recreate invalid analysis fixture");
    fputs("invalid", broken);
    fclose(broken);
    check(ms2_analyze_v1("missing-source.wav", variant) != 0,
        "failed regeneration must report the real source error");
    broken = fopen(variant, "rb");
    check(broken != NULL && fgetc(broken) == 'i',
        "failed regeneration must preserve the existing cache");
    fclose(broken);
    remove(input);
    remove(variant);
    check(exists(profile), "profiling log was not created");
    FILE* log = fopen(profile, "rb");
    char line[8192];
    int hit = 0, miss = 0, skipped = 0, forced = 0, converted = 0;
    while (fgets(line, sizeof(line), log)) {
        hit += strstr(line, "stage=analysis.cache_hit ") != NULL;
        miss += strstr(line, "stage=analysis.missing_begin ") != NULL;
        skipped += strstr(line, "stage=phone.pcm_skipped ") != NULL;
        forced += strstr(line, "stage=analysis.forced_begin ") != NULL;
        converted += strstr(line, "stage=analysis.convert_layer1 ") != NULL;
    }
    fclose(log);
    check(hit == 3 && miss == 1 && skipped == 3, "profiling must distinguish reuse, generation, and skipped synthesis");
    check(forced == 6 && converted == 6,
        "analysis must run for missing/invalid caches and explicit overwrites only");
    remove(profile);
    return 0;
}
