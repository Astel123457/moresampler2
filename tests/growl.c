#define MORESAMPLER2_NO_MAIN
#include "../moresampler2.c"

static void check(int ok, const char* message) {
  if (!ok) { fprintf(stderr, "%s\n", message); exit(1); }
}

static llsm_chunk* model(void) {
  llsm_aoptions* options = llsm_create_aoptions();
  options->thop = .005f;
  llsm_container* conf = llsm_aoptions_toconf(options, 22050);
  llsm_container_attach(conf, LLSM_CONF_NFRM, llsm_create_int(200), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(conf, LLSM_CONF_NSPEC, llsm_create_int(64), llsm_delete_int, llsm_copy_int);
  llsm_chunk* chunk = llsm_create_chunk(conf, 1);
  llsm_delete_container(conf);
  llsm_delete_aoptions(options);
  for (int i = 0; i < 200; ++i) {
    llsm_container_attach(chunk->frames[i], LLSM_FRAME_F0, llsm_create_fp(200), llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(chunk->frames[i], LLSM_FRAME_RD, llsm_create_fp(1), llsm_delete_fp, llsm_copy_fp);
    FP_TYPE* vt = llsm_create_fparray(64);
    for (int j = 0; j < 64; ++j) vt[j] = -10;
    llsm_container_attach(chunk->frames[i], LLSM_FRAME_VTMAGN, vt, llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(chunk->frames[i], LLSM_FRAME_VSPHSE, llsm_create_fparray(64), llsm_delete_fparray, llsm_copy_fparray);
  }
  return chunk;
}

static float* test_signal(int n) {
  float* data = malloc(sizeof(float) * n);
  for (int i = 0; i < n; ++i)
    data[i] = .2f * sin(2 * M_PI * 200 * i / 44100.0)
      + .1f * sin(2 * M_PI * 2400 * i / 44100.0);
  return data;
}

int main(void) {
  llsm_chunk* chunk = model();
  Flags flags;
  parse_flag_string("MG70MD40", &flags);
  ms2_resampler_request request = {0};
  llsm_chunk* flagged = llsm_copy_chunk(chunk);
  apply_remaining_effects(flagged, &request, &flags);
  parse_flag_string("", &flags);
  float growl[] = {70, 70}, distortion[] = {40, 40};
  request.growl_curve = growl; request.distortion_curve = distortion; request.expression_count = 2;
  llsm_chunk* curved = llsm_copy_chunk(chunk);
  apply_remaining_effects(curved, &request, &flags);
  for (int i = 0; i < 200; ++i) {
    for (int key = MS2_FRAME_GROWL; key <= MS2_FRAME_DISTORTION; ++key)
      check(fabs(*(FP_TYPE*)llsm_container_get(flagged->frames[i], key) -
        *(FP_TYPE*)llsm_container_get(curved->frames[i], key)) < 1e-4, "flags and curves differ");
    FP_TYPE* vt = llsm_container_get(curved->frames[i], LLSM_FRAME_VTMAGN);
    check(vt[0] == -10 && vt[63] == -10, "MG/MD must not apply spectral tremolo");
  }
  llsm_delete_chunk(flagged); llsm_delete_chunk(curved);
  llsm_nmframe* nm = llsm_container_get(chunk->frames[0], LLSM_FRAME_NM);
  for (int j = 0; j < nm->npsd; ++j) nm->psd[j] = -100;
  apply_coarseness(chunk, chunk->frames[0], 80);
  int colored = 0;
  for (int j = 0; j < nm->npsd; ++j) {
    check(isfinite(nm->psd[j]), "MC produced nonfinite noise");
    colored |= nm->psd[j] != -100;
  }
  check(colored, "MC must color noise");
  check(!llsm_container_get(chunk->frames[0], MS2_FRAME_GROWL), "MC must not enable MG");
  float* original = test_signal(44100);
  llsm_output output = {0};
  output.fs = 44100; output.ny = 44100; output.y = test_signal(44100);
  check(!ms2_pulse_effects(chunk, &output), "neutral processing failed");
  check(!memcmp(original, output.y, sizeof(float) * 44100), "neutral audio changed");
  float* mg = NULL;
  for (int key = MS2_FRAME_GROWL; key <= MS2_FRAME_DISTORTION; ++key) {
    for (int i = 0; i < 200; ++i) {
      llsm_container_remove(chunk->frames[i], MS2_FRAME_GROWL);
      llsm_container_remove(chunk->frames[i], MS2_FRAME_DISTORTION);
      llsm_container_attach(chunk->frames[i], key, llsm_create_fp(80), llsm_delete_fp, llsm_copy_fp);
    }
    memcpy(output.y, original, sizeof(float) * 44100);
    check(!ms2_pulse_effects(chunk, &output), "active processing failed");
    double difference = 0;
    for (int i = 0; i < 44100; ++i) {
      check(isfinite(output.y[i]) && fabs(output.y[i]) < 2, "unstable pulse output");
      difference += fabs(output.y[i] - original[i]);
    }
    check(difference > 1, "effect is inaudibly small");
    llsm_output short_output = {0};
    short_output.fs = 44100; short_output.ny = 22050; short_output.y = test_signal(22050);
    check(!ms2_pulse_effects(chunk, &short_output), "short processing failed");
    for (int i = 500; i < 21000; ++i)
      check(fabs(short_output.y[i] - output.y[i]) < 1e-6, "effect depends on note length");
    free(short_output.y);
    if (key == MS2_FRAME_GROWL) {
      mg = malloc(sizeof(float) * 44100); memcpy(mg, output.y, sizeof(float) * 44100);
    } else check(memcmp(mg, output.y, sizeof(float) * 44100), "MG and MD must be distinct");
  }
  // The effects must survive the standalone resampler -> wavtool sidecar.
  llsm_aoptions* options = llsm_create_aoptions();
  int fs = 44100, bits = 16, count = 0;
  check(!save_llsm(chunk, "growl-test.llsm2", options, &fs, &bits, 0), "effect cache save failed");
  llsm_chunk* copy = read_llsm("growl-test.llsm2", &count, &fs, &bits, NULL);
  check(copy && count == 200, "effect cache read failed");
  FP_TYPE* value = llsm_container_get(copy->frames[0], MS2_FRAME_DISTORTION);
  check(value && *value == 80, "effect cache lost MD");
  llsm_container_remove(copy->frames[0], MS2_FRAME_DISTORTION);
  ms2_interp_llsm_frame(copy->frames[0], chunk->frames[0], .5f);
  value = llsm_container_get(copy->frames[0], MS2_FRAME_DISTORTION);
  check(value && *value == 40, "stitch did not blend MD");
  llsm_delete_chunk(copy); llsm_delete_aoptions(options);
  remove("growl-test.llsm2");
  free(mg); free(original); free(output.y); llsm_delete_chunk(chunk);
  return 0;
}
