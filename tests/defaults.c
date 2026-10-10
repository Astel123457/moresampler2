#define MORESAMPLER2_NO_MAIN
#include "../moresampler2.c"

static void check(int ok, const char* message) {
  if (!ok) { fprintf(stderr, "%s\n", message); exit(1); }
}

int main(void) {
  Flags flags;
  parse_flag_string("", &flags);
  check(flags.P == 0 && flags.Y == 100, "empty flags must be neutral");
  ms2_resampler_request request = {0};
  check(!remaining_effects_active(&request, &flags), "absent curves must bypass effects");
  float zero[] = {0, 0, 0}, voice[] = {100, 100, 100};
  const float** curves[] = {
    &request.breathiness_curve, &request.mouth_curve, &request.formant_curve,
    &request.singer_formant_curve, &request.growl_curve,
    &request.growl_coarseness_curve, &request.distortion_curve,
    &request.pitch_randomness_curve, &request.unvoiced_emphasis_curve,
    &request.amplitude_modulation_curve, &request.dryness_curve };
  request.expression_count = 3;
  for (size_t i = 0; i < sizeof(curves) / sizeof(curves[0]); ++i) *curves[i] = zero;
  request.voicing_curve = voice;
  request.stability_curve = zero;
  check(!remaining_effects_active(&request, &flags), "flat neutral curves must bypass effects");

  llsm_aoptions* options = llsm_create_aoptions();
  llsm_container* conf = llsm_aoptions_toconf(options, 22050);
  llsm_container_attach(conf, LLSM_CONF_NFRM, llsm_create_int(3), llsm_delete_int, llsm_copy_int);
  llsm_chunk* chunk = llsm_create_chunk(conf, 1);
  llsm_delete_container(conf);
  llsm_delete_aoptions(options);
  FP_TYPE original[] = {-1.1234567f, -20.7654321f, -39.456789f, -10.345678f};
  FP_TYPE* spectra[3];
  for (int i = 0; i < 3; ++i) {
    spectra[i] = llsm_create_fparray(4);
    memcpy(spectra[i], original, sizeof(original));
    llsm_container_attach(chunk->frames[i], LLSM_FRAME_VTMAGN, spectra[i], llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(chunk->frames[i], LLSM_FRAME_RD, llsm_create_fp(1.6f), llsm_delete_fp, llsm_copy_fp);
  }
  apply_gender(chunk, 0);
  apply_gender_curve(chunk, zero, 3);
  apply_tension(chunk, 0);
  apply_tension_curve(chunk, zero, 3);
  apply_remaining_effects(chunk, &request, &flags);
  apply_stabilization(chunk, &request, &flags);
  for (int i = 0; i < 3; ++i) {
    check(llsm_container_get(chunk->frames[i], LLSM_FRAME_VTMAGN) == spectra[i], "neutral effects must not replace spectra");
    check(!memcmp(spectra[i], original, sizeof(original)), "neutral effects must preserve spectra bit-for-bit");
    check(*(FP_TYPE*)llsm_container_get(chunk->frames[i], LLSM_FRAME_RD) == 1.6f,
      "neutral tension must preserve source shape");
  }
  llsm_chunk* tension_curve = llsm_copy_chunk(chunk);
  Flags tension_flags;
  parse_flag_string("Mt100", &tension_flags);
  apply_tension(chunk, tension_flags.Mt);
  float full_tension[] = {100, 100, 100};
  apply_tension_curve(tension_curve, full_tension, 3);
  for (int i = 0; i < 3; ++i) {
    FP_TYPE* rd = llsm_container_get(chunk->frames[i], LLSM_FRAME_RD);
    FP_TYPE* curve_rd = llsm_container_get(tension_curve->frames[i], LLSM_FRAME_RD);
    check(fabs(*rd - .4f) < 1e-6, "full tension must quarter Rd before clamping");
    check(fabs(*rd - *curve_rd) < 1e-6, "Mt and TENC must use the same scaling");
    check(fabs(spectra[i][0] - (original[0] - 6)) < 1e-5 &&
      fabs(spectra[i][3] - (original[3] + 6)) < 1e-5, "tension spectral tilt changed");
    *rd = .4f;
    apply_tension_frame(chunk->frames[i], 100);
    check(*rd == .3f, "positive tension must retain lower Rd limit");
    *rd = 1.6f;
    apply_tension_frame(chunk->frames[i], -100);
    check(*rd == 3.0f, "negative tension must retain upper Rd limit");
  }
  llsm_delete_chunk(tension_curve);
  FP_TYPE waveform[] = {-.123f, .234f, -.345f};
  FP_TYPE before[3]; memcpy(before, waveform, sizeof(before));
  check(normalize_waveform(waveform, 3, .5f, flags.P) == 0, "default normalization gain must be zero");
  check(!memcmp(before, waveform, sizeof(before)), "default normalization must not touch audio");

  float active[] = {0, 25, 0};
  for (size_t i = 0; i < sizeof(curves) / sizeof(curves[0]); ++i) {
    *curves[i] = active;
    check(remaining_effects_active(&request, &flags), "non-neutral curve must enable effects");
    *curves[i] = zero;
  }
  voice[1] = 0;
  check(remaining_effects_active(&request, &flags), "explicit zero voicing is not neutral");
  voice[1] = 100;
  parse_flag_string("Y0", &flags);
  check(remaining_effects_active(&request, &flags), "explicit Y0 must remain active");
  parse_flag_string("P86", &flags);
  check(flags.P == 86, "explicit normalization must remain available");
  llsm_delete_chunk(chunk);
  return 0;
}
