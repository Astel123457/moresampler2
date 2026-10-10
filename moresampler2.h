#ifndef MORESAMPLER2_H
#define MORESAMPLER2_H

#if defined(_WIN32)
#  if defined(MORESAMPLER2_BUILD_SHARED)
#    define MS2_API __declspec(dllexport)
#  else
#    define MS2_API __declspec(dllimport)
#  endif
#else
#  define MS2_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float f0_min, f0_max, voicing_bias;
    int transition_range, refine_f0;
} ms2_analysis_settings;

typedef struct {
    const char* input;
    const char* output;
    float tone;
    float velocity;
    const char* flags;
    float offset;
    float length;
    float consonant;
    float cutoff;
    int volume;
    int modulation;
    int tempo;
    const char* pitch_curve;
    const float* pitch_values;
    int pitch_value_count;
    const float* gender_curve;
    const float* tension_curve;
    const float* breathiness_curve;
    const float* voicing_curve;
    const float* mouth_curve;
    const float* formant_curve;
    const float* singer_formant_curve;
    const float* growl_curve;
    const float* growl_coarseness_curve;
    const float* distortion_curve;
    const float* pitch_randomness_curve;
    const float* unvoiced_emphasis_curve;
    const float* amplitude_modulation_curve;
    const float* stability_curve;
    const float* dryness_curve;
    int expression_count;
    const char* analysis_file; // Optional settings-specific source analysis path.
} ms2_resampler_request;

typedef struct {
    ms2_resampler_request synthesis;
    float position_ms;
    float skip_ms;
    float length_ms;
    float overlap_ms;
    const float* envelope_x_ms;
    const float* envelope_y;
    int envelope_count;
} ms2_phrase_phone;

typedef struct {
    const ms2_phrase_phone* phones;
    int phone_count;
    const char* output;
    float duration_ms;
    int sample_rate;
} ms2_phrase_request;

typedef struct {
    float* samples;
    int sample_count;
    int sample_rate;
} ms2_render_result;

/* Renders one phoneme to request->output. Returns zero on success. */
MS2_API int ms2_resample(const ms2_resampler_request* request);

/* Whole-phrase renderer API. Curves are sampled once per synthesis frame. */
MS2_API int ms2_render_phrase_v1(
    const ms2_phrase_request* request,
    ms2_render_result* result);

/* Generates or replaces reusable analysis for one source audio file. */
MS2_API int ms2_analyze_v1(const char* input, const char* output);
MS2_API int ms2_analyze_config_v1(const char* input, const char* output,
    const ms2_analysis_settings* settings);
MS2_API int ms2_generate_v1(const char* input, const char* output,
    const ms2_analysis_settings* settings, int overwrite);

/* Validates that an analysis file can be read by this engine version. */
MS2_API int ms2_validate_analysis_v1(const char* analysis);

/* Frees engine-owned memory returned by future API additions. */
MS2_API void ms2_free(void* pointer);

// Optional diagnostics; set once before starting worker calls.
MS2_API void ms2_set_profile_path_v1(const char* path);

#ifdef __cplusplus
}
#endif

#endif
