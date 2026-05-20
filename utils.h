#include <stdio.h>
#include <llsm.h>

#define LOG2DB (20.0 / 2.3025851)
#define mag2db(x) (log_2(x) * LOG2DB)

FP_TYPE ms2_lerp(FP_TYPE a, FP_TYPE b, FP_TYPE t);
FP_TYPE ms2_linterpc(FP_TYPE a, FP_TYPE b, FP_TYPE ratio);
void ms2_interp_nmframe(llsm_nmframe* dst, llsm_nmframe* src, FP_TYPE ratio);
void ms2_interp_llsm_frame(llsm_container* dst, llsm_container* src, FP_TYPE ratio);

int write_conf(FILE* f, llsm_aoptions* conf);
int read_conf(FILE* f, llsm_aoptions* opt);

// Saves an llsm_chunk to an .llsm2 file. Returns 0 on success, -1 on failure.
int save_llsm(llsm_chunk* chunk, const char* filename, llsm_aoptions* conf, int* fs, int* nbit);

// Reads an .llsm2 file and returns an llsm_chunk. Returns NULL on failure.
llsm_chunk* read_llsm(const char* filename, int* nfrm, int* fs, int* nbit);

// Gets an llsm_chunk, first reading if an llsm2 file exists, or analyzing the input audio if not. Returns NULL on failure.
llsm_chunk* get_chunk_from_file(const char* filename, int* nfrm, int* fs, int* nbit, int* nx, int* nhop, llsm_aoptions* opt_a_out, llsm_soptions* opt_s_out);