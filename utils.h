#include <stdio.h>
#include <llsm.h>
#include "moresampler2.h"

double ms2_profile_time(void);
void ms2_profile(const char* stage, double started, const char* detail);
int ms2_chunk_is_layer1(llsm_chunk* chunk);
int ms2_prepare_layer1(llsm_chunk* chunk);
llsm_chunk* read_llsm_options(const char* filename, int* nfrm, int* fs, int* nbit,
    float* gain, llsm_aoptions* options);
int ms2_save_analysis(llsm_chunk* chunk, const char* filename, llsm_aoptions* options,
    int* fs, int* bits);

#define LOG2DB (20.0 / 2.3025851)
#define mag2db(x) (log_2(x) * LOG2DB)

FP_TYPE ms2_lerp(FP_TYPE a, FP_TYPE b, FP_TYPE t);
FP_TYPE ms2_linterpc(FP_TYPE a, FP_TYPE b, FP_TYPE ratio);
void ms2_interp_nmframe(llsm_nmframe* dst, llsm_nmframe* src, FP_TYPE ratio);
void ms2_interp_llsm_frame(llsm_container* dst, llsm_container* src, FP_TYPE ratio);

int write_conf(FILE* f, llsm_aoptions* conf);
int read_conf(FILE* f, llsm_aoptions* opt);

// Saves an llsm_chunk to an .llsm2 file. Returns 0 on success, -1 on failure.
int save_llsm(llsm_chunk* chunk, const char* filename, llsm_aoptions* conf, int* fs, int* nbit, float alias_gain_db);

float normalize_waveform(FP_TYPE *waveform, int length, FP_TYPE target_peak, int P_flag);

char* build_llsm_path(const char* filename);

// Reads an .llsm2 file and returns an llsm_chunk. Returns NULL on failure.
llsm_chunk* read_llsm(const char* filename, int* nfrm, int* fs, int* nbit, float* alias_gain_db);

// Gets an llsm_chunk, first reading if an llsm2 file exists, or analyzing the input audio if not. Returns NULL on failure.
llsm_chunk* get_chunk_from_file(const char* filename, int* nfrm, int* fs, int* nbit,
    int* nx, int* nhop, int nosave, int force_analysis, llsm_aoptions* opt_a_out);
llsm_chunk* get_chunk_config(const char* filename, int* nfrm, int* fs, int* nbit,
    int* nx, int* nhop, int nosave, int force_analysis, llsm_aoptions* options,
    const ms2_analysis_settings* settings);
