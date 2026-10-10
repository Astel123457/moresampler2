// Include the implementation to exercise the internal mapper without adding ABI.
#define MORESAMPLER2_NO_MAIN
#include "../moresampler2.c"

static void check(int ok, const char* message) {
  if (!ok) { fprintf(stderr, "%s\n", message); exit(1); }
}

static llsm_chunk* source_chunk(int count) {
  llsm_aoptions* options = llsm_create_aoptions();
  llsm_container* conf = llsm_aoptions_toconf(options, 22050);
  llsm_container_attach(conf, LLSM_CONF_NFRM, llsm_create_int(count),
    llsm_delete_int, llsm_copy_int);
  llsm_chunk* chunk = llsm_create_chunk(conf, 1);
  llsm_delete_container(conf);
  llsm_delete_aoptions(options);
  for (int i = 0; i < count; ++i) {
    llsm_container* frame = chunk->frames[i];
    llsm_container_attach(frame, LLSM_FRAME_F0, llsm_create_fp(100 + i),
      llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(frame, LLSM_FRAME_RD, llsm_create_fp(1),
      llsm_delete_fp, llsm_copy_fp);
    FP_TYPE* spectrum = llsm_create_fparray(4);
    for (int j = 0; j < 4; ++j) spectrum[j] = i;
    llsm_container_attach(frame, LLSM_FRAME_VTMAGN, spectrum,
      llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(frame, LLSM_FRAME_VSPHSE, llsm_create_fparray(1),
      llsm_delete_fparray, llsm_copy_fparray);
  }
  return chunk;
}

static void expect(llsm_chunk* chunk, int frame, double position) {
  FP_TYPE* spectrum = llsm_container_get(chunk->frames[frame], LLSM_FRAME_VTMAGN);
  if (!spectrum || fabs(spectrum[0] - position) > 1e-4) {
    fprintf(stderr, "frame %d: expected source %.6f, got %.6f\n", frame,
      position, spectrum ? spectrum[0] : -999.0); exit(1);
  }
}

int main(void) {
  llsm_chunk* source = source_chunk(101);
  // Reference stretch knots: output [0, 40, 120] -> source [10, 30, 90].
  llsm_chunk* slow = map_resampler_frames(source, 120, 10, 90, 20, 2, 120, 0);
  expect(slow, 1, 10.5); expect(slow, 39, 29.5);
  expect(slow, 40, 30); expect(slow, 80, 60); expect(slow, 119, 89.25);
  llsm_delete_chunk(slow);
  llsm_chunk* neutral = map_resampler_frames(source, 80, 10, 90, 20, 1, 80, 0);
  for (int i = 0; i < 80; ++i) expect(neutral, i, 10 + i);
  llsm_delete_chunk(neutral);
  // Fast attack must still have access to source frames beyond output length.
  llsm_chunk* fast = map_resampler_frames(source, 30, 10, 90, 20, .5, 30, 0);
  expect(fast, 5, 20); expect(fast, 10, 30); expect(fast, 29, 87);
  llsm_delete_chunk(fast);
  // Do not cap the slow consonant to make room for a vowel or lengthen the note.
  llsm_chunk* short_note = map_resampler_frames(source, 10, 10, 90, 20, 2, 10, 0);
  check(*(int*)llsm_container_get(short_note->conf, LLSM_CONF_NFRM) == 10, "length changed");
  expect(short_note, 9, 14.5); llsm_delete_chunk(short_note);
  llsm_chunk* zero = map_resampler_frames(source, 40, 10, 90, 0, 2, 40, 0);
  expect(zero, 0, 10); expect(zero, 39, 88); llsm_delete_chunk(zero);
  llsm_chunk* all = map_resampler_frames(source, 60, 10, 30, 20, .5, 60, 0);
  expect(all, 9, 28); expect(all, 59, 30); llsm_delete_chunk(all);
  // Looping uses alternating source endpoints, as in the original decomp.
  llsm_chunk* loop = map_resampler_frames(source, 200, 10, 90, 20, .5, 200, 1);
  expect(loop, 10, 30); expect(loop, 70, 90); expect(loop, 100, 60);
  expect(loop, 130, 30); expect(loop, 190, 90); llsm_delete_chunk(loop);
  llsm_chunk* fractional = map_resampler_frames(source, 10, 10.25, 90, 2.5, 2, 10, 0);
  expect(fractional, 1, 10.75); expect(fractional, 5, 12.75);
  llsm_delete_chunk(fractional);
  llsm_chunk* last = map_resampler_frames(source, 1, 200, 300, 0, 1, 1, 0);
  expect(last, 0, 100); llsm_delete_chunk(last);
  expect(source, 10, 10); // Mapping must not mutate the source.
  llsm_delete_chunk(source);
  puts("Velocity mapping regressions passed.");
  return 0;
}
