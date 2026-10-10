#include <llsm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libpyin/pyin.h>
#include <ciglet/ciglet.h>
#include <stdint.h>
#include <math.h>
#include <ctype.h>
#include "utils.h"
#include "moresampler2.h"
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#else
#include <time.h>
#include <unistd.h>
#endif

static char ms2_profile_path[4096];

static FP_TYPE pulse_random(uint32_t* state) {
  *state = *state * 1664525u + 1013904223u;
  return (*state >> 8) / (FP_TYPE)16777215;
}

// Adapted from the raw 0040e900/00418b30 trace: pitch-synchronous Hann
// windows, six-pulse MG pattern, alternating MD gain, split at 1500 Hz.
// libllsm2 supplies F0 rather than the original LMIX pulse-position list.
int ms2_pulse_effects(llsm_chunk* chunk, llsm_output* output) {
  if (!chunk || !output || !output->y || output->ny < 2) return 0;
  int* count_ptr = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  FP_TYPE* hop_ptr = llsm_container_get(chunk->conf, LLSM_CONF_THOP);
  if (!count_ptr || *count_ptr <= 0 || !hop_ptr || !isfinite(*hop_ptr) || *hop_ptr <= 0
      || !isfinite(output->fs) || output->fs <= 0) return 1;
  int count = *count_ptr;
  FP_TYPE hop = *hop_ptr;
  int active = 0;
  for (int i = 0; i < count; ++i) {
    FP_TYPE* g = llsm_container_get(chunk->frames[i], MS2_FRAME_GROWL);
    FP_TYPE* d = llsm_container_get(chunk->frames[i], MS2_FRAME_DISTORTION);
    active |= (g && *g > 0) || (d && *d > 0);
  }
  if (!active) return 0; // Neutral synthesis remains bit-for-bit untouched.
  int n = output->ny;
  FP_TYPE* low = calloc((size_t)n, sizeof(FP_TYPE));
  FP_TYPE* result = calloc((size_t)n, sizeof(FP_TYPE));
  FP_TYPE* weights = calloc((size_t)n, sizeof(FP_TYPE));
  if (!low || !result || !weights) {
    free(low); free(result); free(weights); return 1;
  }
  FP_TYPE taps[32], sum = 0;
  for (int k = 0; k < 32; ++k) {
    double x = k - 15.5, fc = 1500.0 / output->fs;
    taps[k] = (FP_TYPE)(sin(2 * M_PI * fc * x) / (M_PI * x)
      * (0.5 - 0.5 * cos(2 * M_PI * k / 31)));
    sum += taps[k];
  }
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < 32; ++k) {
      // Half-sample compensated, centered FIR. High = original - low.
      int j = i + k - 16;
      FP_TYPE a = j >= 0 && j < n ? output->y[j] : 0;
      FP_TYPE b = j + 1 >= 0 && j + 1 < n ? output->y[j + 1] : 0;
      low[i] += (a + b) * (FP_TYPE)0.5 * taps[k] / sum;
    }
  const FP_TYPE pattern[6] = {-.89830506f, -.69491524f, -.1864407f,
    .62711865f, .8305084f, .32203388f};
  uint32_t random = 1;
  double phase = 0;
  unsigned pulse = 0;
  for (int center = 0; center < n; ++center) {
    int frame = min(count - 1, (int)(center / (output->fs * hop)));
    FP_TYPE* f0 = llsm_container_get(chunk->frames[frame], LLSM_FRAME_F0);
    if (!f0 || !isfinite(*f0) || *f0 <= 0) { phase = 0; continue; }
    phase += *f0 / output->fs;
    if (phase < 1) continue;
    phase -= floor(phase);
    FP_TYPE* gp = llsm_container_get(chunk->frames[frame], MS2_FRAME_GROWL);
    FP_TYPE* dp = llsm_container_get(chunk->frames[frame], MS2_FRAME_DISTORTION);
    FP_TYPE g = gp ? *gp / 100 : 0, d = dp ? *dp / 100 : 0;
    int radius = max(1, (int)round(min(n, output->fs / *f0)));
    for (int band = 0; band < 2; ++band) {
      FP_TYPE strength = band ? .5f : .1f;
      int shift = band ? (int)round(.03 * g * radius *
        (cos(pulse * M_PI / 3) - .5 + pulse_random(&random))) : 0;
      FP_TYPE gain = 1 + 1.5f * g * ((pulse_random(&random) - .5f) * .3f
        + strength * pattern[pulse % 6]);
      gain += strength * 1.5f * d * ((pulse & 1 ? -1.5f : .5f) + pulse_random(&random));
      for (int j = -radius; j < radius; ++j) {
        int src = center + j, dst = src + shift;
        if (src < 0 || src >= n) continue;
        FP_TYPE w = (FP_TYPE)(.5 - .5 * cos(2 * M_PI * (j + radius) / (2 * radius - 1)));
        if (!band) weights[src] += w;
        if (dst >= 0 && dst < n)
          result[dst] += w * gain * (band ? output->y[src] - low[src] : low[src]);
      }
    }
    ++pulse;
  }
  for (int i = 0; i < n; ++i)
    output->y[i] = result[i] + max(0, 1 - weights[i]) * output->y[i];
  free(low); free(result); free(weights);
  return 0;
}

// Configure once, before rendering/validation starts. Never changes audio or caches.
void ms2_set_profile_path_v1(const char* path) {
  if (!ms2_profile_path[0] && path && strlen(path) < sizeof(ms2_profile_path))
    strcpy(ms2_profile_path, path);
}

double ms2_profile_time(void) {
#ifdef _WIN32
  LARGE_INTEGER now, frequency;
  QueryPerformanceCounter(&now);
  QueryPerformanceFrequency(&frequency);
  return (double)now.QuadPart / frequency.QuadPart;
#else
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return now.tv_sec + now.tv_nsec * 1e-9;
#endif
}

void ms2_profile(const char* stage, double started, const char* detail) {
  if (!ms2_profile_path[0]) return;
  double now = ms2_profile_time();
  FILE* file = fopen(ms2_profile_path, "a");
  if (!file) return; // Profiling must never make rendering fail.
  fprintf(file, "native t=%.6f stage=%s elapsed_ms=%.3f detail=%s\n",
    now, stage, (now - started) * 1000, detail ? detail : "");
  fclose(file);
}

#define MS2_CACHE_VERSION_RAW 2
#define MS2_CACHE_VERSION_EFFECTS 4
#define MS2_CACHE_VERSION_L1_HYBRID 3
#define MS2_VTMAGN_COEFFS 384
#define MS2_NPSD_COEFFS 96

static void ddct(int n, int isgn, FP_TYPE *a);

static FP_TYPE* ms2_encode_dct(FP_TYPE* values, int n, int order) {
  int coeff_count = min(order, n);
  FP_TYPE* dct_full = dct(values, n);
  FP_TYPE* coeffs = llsm_create_fparray(coeff_count);
  if (!dct_full || !coeffs) {
    free(dct_full);
    llsm_delete_fparray(coeffs);
    return NULL;
  }
  memcpy(coeffs, dct_full, sizeof(FP_TYPE) * (size_t)coeff_count);
  free(dct_full);
  return coeffs;
}

static FP_TYPE* ms2_decode_dct(FP_TYPE* coeffs, int coeff_count, int n) {
  FP_TYPE* values = llsm_create_fparray(n);
  if (!values) return NULL;

  FP_TYPE inv_sqrt_n = (FP_TYPE)(1.0 / sqrt((double)n));
  FP_TYPE inv_sqrt_half_n = (FP_TYPE)sqrt(2.0 / (double)n);
  for (int i = 0; i < n; ++i) {
    FP_TYPE value = coeffs[0] * inv_sqrt_n;
    for (int k = 1; k < coeff_count; ++k) {
      value += coeffs[k] * inv_sqrt_half_n *
               (FP_TYPE)cos((double)M_PI / (double)n * (double)(i + 0.5) * (double)k);
    }
    values[i] = value;
  }
  return values;
}

float normalize_waveform(FP_TYPE *waveform, int length, FP_TYPE target_peak,
                        int P_flag) {
  if (P_flag <= 0)
    return 0.0f;

  FP_TYPE peak = 0.0f;
  for (int i = 0; i < length; ++i) {
    FP_TYPE abs_val = fabsf(waveform[i]);
    if (abs_val > peak)
      peak = abs_val;
  }

  if (peak < 1e-9f)
    return 0.0f;

  FP_TYPE full_scale = target_peak / peak;
  FP_TYPE blend = P_flag / 100.0f;
  FP_TYPE scale = linterp(1.0f, full_scale, blend);

  for (int i = 0; i < length; ++i)
    waveform[i] *= scale;
  return 20.0f * log10f(scale);
}

FP_TYPE ms2_lerp(FP_TYPE a, FP_TYPE b, FP_TYPE t) {
  return a + (b - a) * t;
}

FP_TYPE ms2_linterpc(FP_TYPE a, FP_TYPE b, FP_TYPE ratio) {
  FP_TYPE ax = (FP_TYPE)cos((double)a);
  FP_TYPE ay = (FP_TYPE)sin((double)a);
  FP_TYPE bx = (FP_TYPE)cos((double)b);
  FP_TYPE by = (FP_TYPE)sin((double)b);
  FP_TYPE cx = ms2_lerp(ax, bx, ratio);
  FP_TYPE cy = ms2_lerp(ay, by, ratio);
  return (FP_TYPE)atan2((double)cy, (double)cx);
}

void ms2_interp_nmframe(llsm_nmframe* dst, llsm_nmframe* src, FP_TYPE ratio) {
  if (!dst || !src) return;


  int npsd = min(dst->npsd, src->npsd);
  for (int i = 0; i < npsd; ++i) {
    dst->psd[i] = ms2_lerp(dst->psd[i], src->psd[i], ratio);
  }

  int nchan = min(dst->nchannel, src->nchannel);
  for (int b = 0; b < nchan; ++b) {
    llsm_hmframe* srceenv = src->eenv[b];
    llsm_hmframe* dsteenv = dst->eenv[b];
    if (!srceenv || !dsteenv) continue;

    dst->edc[b] = ms2_lerp(dst->edc[b], src->edc[b], ratio);
    int b_minnhar = min(srceenv->nhar, dsteenv->nhar);
    int b_maxnhar = max(srceenv->nhar, dsteenv->nhar);
    if (dsteenv->nhar < b_maxnhar) {
      dsteenv->ampl = (FP_TYPE*)realloc(dsteenv->ampl, sizeof(FP_TYPE) * (size_t)b_maxnhar);
      dsteenv->phse = (FP_TYPE*)realloc(dsteenv->phse, sizeof(FP_TYPE) * (size_t)b_maxnhar);
      if (!dsteenv->ampl || !dsteenv->phse) continue;
    }
    for (int i = 0; i < b_minnhar; ++i) {
      dsteenv->ampl[i] = ms2_lerp(dsteenv->ampl[i], srceenv->ampl[i], ratio);
      dsteenv->phse[i] = ms2_linterpc(dsteenv->phse[i], srceenv->phse[i], ratio);
    }
    if (b_maxnhar == srceenv->nhar) {
      for (int i = b_minnhar; i < b_maxnhar; ++i) {
        dsteenv->ampl[i] = srceenv->ampl[i];
        dsteenv->phse[i] = srceenv->phse[i];
      }
    }
    dsteenv->nhar = b_maxnhar;
  }
}

void ms2_interp_llsm_frame(llsm_container* dst, llsm_container* src, FP_TYPE ratio) {
#define MS2_EPS ((FP_TYPE)1e-8)
#define MS2_LOG2DB ((FP_TYPE)(20.0 / 2.3025851))
  if (!dst || !src) return;
  for (int key = MS2_FRAME_GROWL; key <= MS2_FRAME_DISTORTION; ++key) {
    FP_TYPE* a = llsm_container_get(dst, key);
    FP_TYPE* b = llsm_container_get(src, key);
    if (a || b) llsm_container_attach(dst, key,
      llsm_create_fp((a ? *a : 0) * (1 - ratio) + (b ? *b : 0) * ratio),
      llsm_delete_fp, llsm_copy_fp);
  }

  FP_TYPE* dst_f0_p = (FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_F0);
  FP_TYPE* src_f0_p = (FP_TYPE*)llsm_container_get(src, LLSM_FRAME_F0);
  if (!dst_f0_p || !src_f0_p) return;
  FP_TYPE dst_f0 = *dst_f0_p;
  FP_TYPE src_f0 = *src_f0_p;

  llsm_nmframe* dst_nm = (llsm_nmframe*)llsm_container_get(dst, LLSM_FRAME_NM);
  llsm_nmframe* src_nm = (llsm_nmframe*)llsm_container_get(src, LLSM_FRAME_NM);
  FP_TYPE* src_rd = (FP_TYPE*)llsm_container_get(src, LLSM_FRAME_RD);
  FP_TYPE* dst_rd = (FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_RD);
  FP_TYPE* dst_vsphse = (FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_VSPHSE);
  FP_TYPE* src_vsphse = (FP_TYPE*)llsm_container_get(src, LLSM_FRAME_VSPHSE);
  FP_TYPE* dst_vtmagn = (FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_VTMAGN);
  FP_TYPE* src_vtmagn = (FP_TYPE*)llsm_container_get(src, LLSM_FRAME_VTMAGN);

  llsm_container* voiced = (dst_f0 <= 0 && src_f0 <= 0) ? NULL : (src_f0 > 0 ? src : dst);
  int bothvoiced = (dst_f0 > 0 && src_f0 > 0);

  int dstnhar = dst_vsphse == NULL ? 0 : llsm_fparray_length(dst_vsphse);
  int srcnhar = src_vsphse == NULL ? 0 : llsm_fparray_length(src_vsphse);
  int maxnhar = max(dstnhar, srcnhar);
  int minnhar = min(dstnhar, srcnhar);

  if (!bothvoiced && voiced == src && src_rd) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(src_f0), llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(*src_rd), llsm_delete_fp, llsm_copy_fp);
  } else if (voiced == NULL) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp((FP_TYPE)0), llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp((FP_TYPE)1.0), llsm_delete_fp, llsm_copy_fp);
  }

  int nspec = dst_vtmagn != NULL ? llsm_fparray_length(dst_vtmagn)
                                 : (src_vtmagn != NULL ? llsm_fparray_length(src_vtmagn) : 0);

  if (bothvoiced && dst_rd && src_rd) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(ms2_lerp(dst_f0, src_f0, ratio)), llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(ms2_lerp(*dst_rd, *src_rd, ratio)), llsm_delete_fp, llsm_copy_fp);

    FP_TYPE* vsphse = llsm_create_fparray(maxnhar);
    FP_TYPE* vtmagn = llsm_create_fparray(nspec);
    for (int i = 0; i < minnhar; ++i) {
      vsphse[i] = ms2_linterpc(dst_vsphse[i], src_vsphse[i], ratio);
    }
    for (int i = 0; i < nspec; ++i) {
      vtmagn[i] = ms2_lerp(dst_vtmagn[i], src_vtmagn[i], ratio);
    }
    if (dstnhar < srcnhar) {
      for (int i = minnhar; i < maxnhar; ++i) {
        vsphse[i] = src_vsphse[i];
      }
    }

    llsm_container_attach(dst, LLSM_FRAME_VSPHSE, vsphse, llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(dst, LLSM_FRAME_VTMAGN, vtmagn, llsm_delete_fparray, llsm_copy_fparray);
    dst_vtmagn = (FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_VTMAGN);
  } else if (voiced == src && src_vsphse && src_vtmagn) {
    llsm_container_attach(dst, LLSM_FRAME_VSPHSE, llsm_copy_fparray(src_vsphse), llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(dst, LLSM_FRAME_VTMAGN, llsm_copy_fparray(src_vtmagn), llsm_delete_fparray, llsm_copy_fparray);
    dst_vtmagn = (FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_VTMAGN);
    if (dst_vtmagn) {
      FP_TYPE fade = (FP_TYPE)(log((double)max(MS2_EPS, ratio)) * (double)MS2_LOG2DB);
      for (int i = 0; i < nspec; ++i) dst_vtmagn[i] += fade;
    }
  } else if (dst_vtmagn) {
    FP_TYPE fade = (FP_TYPE)(log((double)max(MS2_EPS, (FP_TYPE)1.0 - ratio)) * (double)MS2_LOG2DB);
    for (int i = 0; i < nspec; ++i) dst_vtmagn[i] += fade;
  }

  if (dst_vtmagn) {
    for (int i = 0; i < nspec; ++i) dst_vtmagn[i] = max((FP_TYPE)-80, dst_vtmagn[i]);
  }

  if (dst_nm && src_nm) {
    ms2_interp_nmframe(dst_nm, src_nm, ratio);
  }
#undef MS2_EPS
#undef MS2_LOG2DB
}

int write_conf(FILE* f, llsm_aoptions* conf) {
    fwrite(&conf->thop, sizeof(FP_TYPE), 1, f);
    fwrite(&conf->maxnhar, sizeof(int), 1, f);
    fwrite(&conf->maxnhar_e, sizeof(int), 1, f);
    fwrite(&conf->npsd, sizeof(int), 1, f);
    fwrite(&conf->nchannel, sizeof(int), 1, f);
    int chanfreq_len = conf->nchannel;
    fwrite(&chanfreq_len, sizeof(int), 1, f);
    fwrite(conf->chanfreq, sizeof(FP_TYPE), chanfreq_len, f);
    fwrite(&conf->lip_radius, sizeof(FP_TYPE), 1, f);
    fwrite(&conf->f0_refine, sizeof(FP_TYPE), 1, f);
    fwrite(&conf->hm_method, sizeof(int), 1, f);
    fwrite(&conf->rel_winsize, sizeof(FP_TYPE), 1, f);
    return 0;
}

int read_conf(FILE* f, llsm_aoptions* opt) {
    fread(&opt->thop, sizeof(FP_TYPE), 1, f);
    fread(&opt->maxnhar, sizeof(int), 1, f);
    fread(&opt->maxnhar_e, sizeof(int), 1, f);
    fread(&opt->npsd, sizeof(int), 1, f);
    fread(&opt->nchannel, sizeof(int), 1, f);
    int chanfreq_len;
    fread(&chanfreq_len, sizeof(int), 1, f);
    free(opt->chanfreq);
    opt->chanfreq = malloc(sizeof(FP_TYPE) * chanfreq_len);
    if (!opt->chanfreq) return -1;
    fread(opt->chanfreq, sizeof(FP_TYPE), chanfreq_len, f);
    fread(&opt->lip_radius, sizeof(FP_TYPE), 1, f);
    fread(&opt->f0_refine, sizeof(FP_TYPE), 1, f);
    fread(&opt->hm_method, sizeof(int), 1, f);
    fread(&opt->rel_winsize, sizeof(FP_TYPE), 1, f);
    return 0;
}

int save_llsm(llsm_chunk* chunk, const char* filename, llsm_aoptions* conf, int* fs, int* nbit, float alias_gain_db) {
  int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  int nfrm_val = nfrm_p ? *nfrm_p : 0;
  int layer = (nfrm_val > 0 && llsm_frame_checklayer1(chunk->frames[0])) ? 1 : 0;
  int nspec = 0;
  int vtmagn_coeffs = MS2_VTMAGN_COEFFS;
  int npsd_coeffs = MS2_NPSD_COEFFS;

  if (layer == 1) {
    int* nspec_p = llsm_container_get(chunk->conf, LLSM_CONF_NSPEC);
    nspec = nspec_p ? *nspec_p : 0;
    vtmagn_coeffs = min(vtmagn_coeffs, nspec);
  }
  npsd_coeffs = min(npsd_coeffs, conf->npsd);

  FILE* f = fopen(filename, "wb");
  if (!f) return -1;

  // Wavtool reloads and resaves layer-1 caches across note concatenation, so
  // the on-disk cache must stay lossless on the write path.
  int version = MS2_CACHE_VERSION_RAW;
  for (int i = 0; i < nfrm_val; ++i)
    if (llsm_container_get(chunk->frames[i], MS2_FRAME_GROWL) ||
        llsm_container_get(chunk->frames[i], MS2_FRAME_DISTORTION))
      version = MS2_CACHE_VERSION_EFFECTS;
  fwrite("LLSM2", 1, 5, f);
  fwrite(&version, sizeof(int), 1, f);
  fwrite(&layer, sizeof(int), 1, f);
  fwrite(&nspec, sizeof(int), 1, f);
  fwrite(&alias_gain_db, sizeof(float), 1, f);
  fwrite(&nfrm_val, sizeof(int), 1, f);
  fwrite(fs, sizeof(int), 1, f);
  fwrite(nbit, sizeof(int), 1, f);
  write_conf(f, conf);
  if (version == MS2_CACHE_VERSION_L1_HYBRID) {
    fwrite(&vtmagn_coeffs, sizeof(int), 1, f);
    fwrite(&npsd_coeffs, sizeof(int), 1, f);
  }

  for (int i = 0; i < nfrm_val; ++i) {
    llsm_container* frame = chunk->frames[i];
    FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
    fwrite(f0, sizeof(FP_TYPE), 1, f);

    if (layer == 0) {
      llsm_hmframe* hm = llsm_container_get(frame, LLSM_FRAME_HM);
      fwrite(&hm->nhar, sizeof(int), 1, f);
      fwrite(hm->ampl, sizeof(FP_TYPE), hm->nhar, f);
      fwrite(hm->phse, sizeof(FP_TYPE), hm->nhar, f);
    } else {
      FP_TYPE* rd = llsm_container_get(frame, LLSM_FRAME_RD);
      FP_TYPE rd_val = rd ? *rd : (FP_TYPE)1.0;
      FP_TYPE* vtmagn = llsm_container_get(frame, LLSM_FRAME_VTMAGN);
      FP_TYPE* vsphse = llsm_container_get(frame, LLSM_FRAME_VSPHSE);
      int vsphse_n = (vsphse && *f0 > 0) ? llsm_fparray_length(vsphse) : 0;
      fwrite(&rd_val, sizeof(FP_TYPE), 1, f);

      if (version == MS2_CACHE_VERSION_L1_HYBRID) {
        if (*f0 > 0) {
          FP_TYPE* vtmagn_coeff = ms2_encode_dct(vtmagn, llsm_fparray_length(vtmagn), vtmagn_coeffs);
          if (!vtmagn_coeff) {
            fclose(f);
            return -1;
          }
          fwrite(vtmagn_coeff, sizeof(FP_TYPE), vtmagn_coeffs, f);
          llsm_delete_fparray(vtmagn_coeff);
        }
      } else {
        int vtmagn_n = (vtmagn && *f0 > 0) ? llsm_fparray_length(vtmagn) : 0;
        fwrite(&vtmagn_n, sizeof(int), 1, f);
        if (vtmagn_n > 0) fwrite(vtmagn, sizeof(FP_TYPE), vtmagn_n, f);
      }

      fwrite(&vsphse_n, sizeof(int), 1, f);
      if (vsphse_n > 0) fwrite(vsphse, sizeof(FP_TYPE), vsphse_n, f);
    }

    llsm_nmframe* nm = llsm_container_get(frame, LLSM_FRAME_NM);
    if (version == MS2_CACHE_VERSION_L1_HYBRID && layer == 1) {
      FP_TYPE* psd_coeff = ms2_encode_dct(nm->psd, nm->npsd, npsd_coeffs);
      if (!psd_coeff) {
        fclose(f);
        return -1;
      }
      fwrite(psd_coeff, sizeof(FP_TYPE), npsd_coeffs, f);
      llsm_delete_fparray(psd_coeff);
    } else {
      fwrite(&nm->npsd, sizeof(int), 1, f);
      fwrite(nm->psd, sizeof(FP_TYPE), nm->npsd, f);
      fwrite(&nm->nchannel, sizeof(int), 1, f);
    }

    for (int j = 0; j < nm->nchannel; ++j) {
      llsm_hmframe* eenv = nm->eenv[j];
      fwrite(&nm->edc[j], sizeof(FP_TYPE), 1, f);
      fwrite(&eenv->nhar, sizeof(int), 1, f);
      fwrite(eenv->ampl, sizeof(FP_TYPE), eenv->nhar, f);
      fwrite(eenv->phse, sizeof(FP_TYPE), eenv->nhar, f);
    }
    if (version == MS2_CACHE_VERSION_EFFECTS) {
      for (int key = MS2_FRAME_GROWL; key <= MS2_FRAME_DISTORTION; ++key) {
        FP_TYPE* value = llsm_container_get(frame, key);
        FP_TYPE scalar = value ? *value : 0;
        fwrite(&scalar, sizeof(scalar), 1, f);
      }
    }
  }

  int failed = ferror(f);
  if (fclose(f) != 0) failed = 1;
  return failed ? -1 : 0;
}

llsm_chunk* read_llsm_options(const char* filename, int* nfrm, int* fs, int* nbit,
    float* alias_gain_db, llsm_aoptions* options) {
  FILE* f = fopen(filename, "rb");
  if (!f) return NULL;

  char header[5];
  fread(header, 1, 5, f);
  if (strncmp(header, "LLSM2", 5) != 0) { fclose(f); return NULL; }

  int version;
  fread(&version, sizeof(int), 1, f);
  if (version != 1 && version != MS2_CACHE_VERSION_RAW && version != MS2_CACHE_VERSION_L1_HYBRID
      && version != MS2_CACHE_VERSION_EFFECTS) {
    fclose(f);
    return NULL;
  }

  int layer = 0;
  int nspec = 0;
  float tmp_alias_gain_db = 0.0f;

  if (version >= 2) {
      fread(&layer, sizeof(int), 1, f);
      fread(&nspec, sizeof(int), 1, f);

      if (alias_gain_db)
          fread(alias_gain_db, sizeof(float), 1, f);
      else
          fread(&tmp_alias_gain_db, sizeof(float), 1, f);
  }

  fread(nfrm, sizeof(int), 1, f);
  fread(fs, sizeof(int), 1, f);
  fread(nbit, sizeof(int), 1, f);

  llsm_aoptions* aopt = options ? options : llsm_create_aoptions();
  if (read_conf(f, aopt) != 0) {
    if (!options) llsm_delete_aoptions(aopt);
    fclose(f);
    return NULL;
  }

  llsm_container* conf = llsm_aoptions_toconf(aopt, *fs / 2.0);
  int conf_npsd = aopt->npsd;
  int conf_nchannel = aopt->nchannel;
  if (!options) llsm_delete_aoptions(aopt);
  llsm_container_attach(conf, LLSM_CONF_NFRM,
    llsm_create_int(*nfrm), llsm_delete_int, llsm_copy_int);
  if (layer == 1 && nspec > 0) {
    llsm_container_attach(conf, LLSM_CONF_NSPEC,
      llsm_create_int(nspec), llsm_delete_int, llsm_copy_int);
  }

  int vtmagn_coeffs = 0;
  int npsd_coeffs = 0;
  if (version == MS2_CACHE_VERSION_L1_HYBRID) {
    fread(&vtmagn_coeffs, sizeof(int), 1, f);
    fread(&npsd_coeffs, sizeof(int), 1, f);
  }

  // Frames are decoded below; the second argument is an initialization flag.
  llsm_chunk* chunk = llsm_create_chunk(conf, 0);
  llsm_delete_container(conf);
  if (!chunk) { fclose(f); return NULL; }

  for (int i = 0; i < *nfrm; ++i) {
    llsm_container* frame = llsm_create_frame(0, 0, 0, 0);
    // Own the frame immediately so decode-error cleanup also releases it.
    chunk->frames[i] = frame;
    FP_TYPE* f0 = malloc(sizeof(FP_TYPE));
    fread(f0, sizeof(FP_TYPE), 1, f);
    llsm_container_attach(frame, LLSM_FRAME_F0, f0, free, llsm_copy_fp);

    if (layer == 0) {
      int nhar;
      fread(&nhar, sizeof(int), 1, f);
      llsm_hmframe* hm = llsm_create_hmframe(nhar);
      fread(hm->ampl, sizeof(FP_TYPE), nhar, f);
      fread(hm->phse, sizeof(FP_TYPE), nhar, f);
      llsm_container_attach(frame, LLSM_FRAME_HM, hm, llsm_delete_hmframe, llsm_copy_hmframe);
    } else {
      llsm_container_remove(frame, LLSM_FRAME_HM);
      FP_TYPE* rd = malloc(sizeof(FP_TYPE));
      fread(rd, sizeof(FP_TYPE), 1, f);
      llsm_container_attach(frame, LLSM_FRAME_RD, rd, llsm_delete_fp, llsm_copy_fp);

      if (version == MS2_CACHE_VERSION_L1_HYBRID && *f0 > 0) {
        FP_TYPE* coeff = llsm_create_fparray(vtmagn_coeffs);
        FP_TYPE* vtmagn = NULL;
        fread(coeff, sizeof(FP_TYPE), vtmagn_coeffs, f);
        vtmagn = ms2_decode_dct(coeff, vtmagn_coeffs, nspec);
        llsm_delete_fparray(coeff);
        if (!vtmagn) {
          llsm_delete_chunk(chunk);
          fclose(f);
          return NULL;
        }
        llsm_container_attach(frame, LLSM_FRAME_VTMAGN, vtmagn, llsm_delete_fparray, llsm_copy_fparray);
      } else {
        int vtmagn_n;
        fread(&vtmagn_n, sizeof(int), 1, f);
        if (vtmagn_n > 0) {
          FP_TYPE* vtmagn = llsm_create_fparray(vtmagn_n);
          fread(vtmagn, sizeof(FP_TYPE), vtmagn_n, f);
          llsm_container_attach(frame, LLSM_FRAME_VTMAGN, vtmagn, llsm_delete_fparray, llsm_copy_fparray);
        }
      }

      int vsphse_n;
      fread(&vsphse_n, sizeof(int), 1, f);
      if (vsphse_n > 0) {
        FP_TYPE* vsphse = llsm_create_fparray(vsphse_n);
        fread(vsphse, sizeof(FP_TYPE), vsphse_n, f);
        llsm_container_attach(frame, LLSM_FRAME_VSPHSE, vsphse, llsm_delete_fparray, llsm_copy_fparray);
      }
    }

    llsm_nmframe* nm = malloc(sizeof(llsm_nmframe));
    if (version == MS2_CACHE_VERSION_L1_HYBRID && layer == 1) {
      FP_TYPE* coeff = llsm_create_fparray(npsd_coeffs);
      nm->npsd = conf_npsd;
      nm->nchannel = conf_nchannel;
      fread(coeff, sizeof(FP_TYPE), npsd_coeffs, f);
      nm->psd = ms2_decode_dct(coeff, npsd_coeffs, nm->npsd);
      llsm_delete_fparray(coeff);
      if (!nm->psd) {
        free(nm);
        llsm_delete_chunk(chunk);
        fclose(f);
        return NULL;
      }
    } else {
      fread(&nm->npsd, sizeof(int), 1, f);
      nm->psd = malloc(sizeof(FP_TYPE) * (size_t)nm->npsd);
      fread(nm->psd, sizeof(FP_TYPE), nm->npsd, f);
      fread(&nm->nchannel, sizeof(int), 1, f);
    }

    nm->edc = malloc(sizeof(FP_TYPE) * (size_t)nm->nchannel);
    nm->eenv = malloc(sizeof(llsm_hmframe*) * (size_t)nm->nchannel);
    for (int j = 0; j < nm->nchannel; ++j) {
      fread(&nm->edc[j], sizeof(FP_TYPE), 1, f);
      int nhar_e;
      fread(&nhar_e, sizeof(int), 1, f);
      llsm_hmframe* eenv = llsm_create_hmframe(nhar_e);
      fread(eenv->ampl, sizeof(FP_TYPE), nhar_e, f);
      fread(eenv->phse, sizeof(FP_TYPE), nhar_e, f);
      nm->eenv[j] = eenv;
    }
    llsm_container_attach(frame, LLSM_FRAME_NM, nm, llsm_delete_nmframe, llsm_copy_nmframe);
    if (version == MS2_CACHE_VERSION_EFFECTS) {
      for (int key = MS2_FRAME_GROWL; key <= MS2_FRAME_DISTORTION; ++key) {
        FP_TYPE value;
        if (fread(&value, sizeof(value), 1, f) != 1 || !isfinite(value) || value < 0 || value > 100) {
          llsm_delete_chunk(chunk); fclose(f); return NULL;
        }
        if (value > 0) llsm_container_attach(frame, key,
          llsm_create_fp(value), llsm_delete_fp, llsm_copy_fp);
      }
    }
  }

  fclose(f);
  return chunk;
}

llsm_chunk* read_llsm(const char* filename, int* nfrm, int* fs, int* nbit, float* gain) {
  return read_llsm_options(filename, nfrm, fs, nbit, gain, NULL);
}

int ms2_chunk_is_layer1(llsm_chunk* chunk) {
  if (!chunk) return 0;
  int* nspec = llsm_container_get(chunk->conf, LLSM_CONF_NSPEC);
  if (!nspec || *nspec <= 0) return 0;
  int* count = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  if (!count || *count <= 0) return 0;
  for (int i = 0; i < *count; ++i)
    if (!chunk->frames[i] || !llsm_frame_checklayer1(chunk->frames[i])) return 0;
  return 1;
}

int ms2_prepare_layer1(llsm_chunk* chunk) {
  if (ms2_chunk_is_layer1(chunk)) return 0;
  double started = ms2_profile_time();
  llsm_chunk_tolayer1(chunk, 2048);
  ms2_profile("analysis.convert_layer1", started, NULL);
  return ms2_chunk_is_layer1(chunk) ? 0 : 1;
}

// Replace only after a complete write, preserving the old cache on failure.
int ms2_save_analysis(llsm_chunk* chunk, const char* filename, llsm_aoptions* options,
    int* fs, int* bits) {
  char* temporary = malloc(strlen(filename) + 8);
  if (!temporary) return -1;
  sprintf(temporary, "%s.XXXXXX", filename);
#ifdef _WIN32
  int fd = _mktemp_s(temporary, strlen(temporary) + 1) == 0
    ? _open(temporary, _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY, _S_IREAD | _S_IWRITE) : -1;
#else
  int fd = mkstemp(temporary);
#endif
  if (fd < 0) { free(temporary); return -1; }
#ifdef _WIN32
  _close(fd);
#else
  close(fd);
#endif
  int result = save_llsm(chunk, temporary, options, fs, bits, 0);
  if (!result) {
#ifdef _WIN32
    result = MoveFileExA(temporary, filename, MOVEFILE_REPLACE_EXISTING) ? 0 : -1;
#else
    result = rename(temporary, filename);
#endif
  }
  if (result) remove(temporary);
  free(temporary);
  return result;
}

char* build_llsm_path(const char* filename) {
    size_t llsm_path_len = strlen(filename) + 7; // ".llsm2" + null terminator
    char* llsm_path = malloc(llsm_path_len);
    if (!llsm_path) return NULL;
    snprintf(llsm_path, llsm_path_len, "%s", filename);
    char* ext = strrchr(llsm_path, '.');
    if (ext) strcpy(ext, ".llsm2"); // Replace extension with .llsm2
    else strcat(llsm_path, ".llsm2"); // No extension, just append
    return llsm_path;
}

llsm_chunk* get_chunk_from_file(const char* filename, int* nfrm, int* fs, int* nbit,
    int* nx, int* nhop, int nosave, int force_analysis, llsm_aoptions* opt_a) {
    return get_chunk_config(filename, nfrm, fs, nbit, nx, nhop, nosave,
        force_analysis, opt_a, NULL);
}

llsm_chunk* get_chunk_config(const char* filename, int* nfrm, int* fs, int* nbit,
    int* nx, int* nhop, int nosave, int force_analysis, llsm_aoptions* opt_a,
    const ms2_analysis_settings* settings) {

    double profile_started = ms2_profile_time();
    llsm_chunk* chunk = NULL;
    FP_TYPE* f0 = NULL;
    FP_TYPE* input = NULL;

    // Build expected .llsm2 path from input WAV path
    char* llsm_path = build_llsm_path(filename);
    if (!llsm_path) {
      printf("WARNING: Failed to build LLSM2 filename path. Analysis will proceed without caching.\n");
    }
    FILE* llsm_file = llsm_path ? fopen(llsm_path, "rb") : NULL;

    if (llsm_file && !force_analysis) {
      // File exists — use cached analysis
      fclose(llsm_file);
      printf("Loading cached LLSM analysis: %s\n", llsm_path);
      chunk = read_llsm_options(llsm_path, nfrm, fs, nbit, NULL, opt_a);
      ms2_profile(chunk ? "analysis.cache_hit" : "analysis.cache_invalid", profile_started, llsm_path);
      if (chunk && !ms2_chunk_is_layer1(chunk)) {
        if (ms2_prepare_layer1(chunk)) {
          llsm_delete_chunk(chunk);
          chunk = NULL;
        }
      }
      free(llsm_path);
      if (!chunk) {
        printf("Failed to read .llsm2 file\n");
        return NULL;
      }
    } else {
      if (llsm_file) fclose(llsm_file);
      // No cache — analyze audio
      ms2_profile(force_analysis ? "analysis.forced_begin" : "analysis.missing_begin", profile_started, filename);
      printf("Reading input WAV: %s\n", filename);
      input = wavread((char*)filename, fs, nbit, nx);
      if (!input) { free(llsm_path); return NULL; }
    
      printf("Estimating F0\n");

      float* input_temp;
      //copy input to input_temp
      input_temp = malloc(sizeof(float) * (size_t)(*nx));
      if (!input_temp) { free(input); free(llsm_path); return NULL; }
      for (int i = 0; i < *nx; ++i) {
        input_temp[i] = input[i];
      }
      normalize_waveform(input_temp, *nx, 0.99f, 100);
      pyin_config param = pyin_init(*nhop);
      param.fmin = settings ? settings->f0_min : 50.0f;
      param.fmax = settings ? settings->f0_max : 800.0f;
      param.trange = settings ? settings->transition_range : 24;
      param.bias = settings ? settings->voicing_bias : 2;
      param.nf = ceil(*fs * 0.025);
      f0 = pyin_analyze(param, input_temp, *nx, *fs, nfrm);
      if (!f0) { free(input); free(input_temp); free(llsm_path); return NULL; }
    
      opt_a->thop = (FP_TYPE)(*nhop) / (*fs);
      opt_a->f0_refine = settings ? settings->refine_f0 : 1;
      opt_a->hm_method = LLSM_AOPTION_HMCZT;
    
      printf("Analysis\n");
      chunk = llsm_analyze(opt_a, input, *nx, *fs, f0, *nfrm, NULL);
      if (chunk && ms2_prepare_layer1(chunk)) { llsm_delete_chunk(chunk); chunk = NULL; }
      if (!chunk) { free(input); free(f0); free(input_temp); free(llsm_path); return NULL; }
      if(llsm_path != NULL) {
        printf("Saving analysis result to cache: %s\n", llsm_path);
          if (!nosave) {
            if (ms2_save_analysis(chunk, llsm_path, opt_a, fs, nbit) != 0) {
              printf("Failed to save .llsm2 file.\n");
            }
          }
      }
      free(llsm_path);
      free(input_temp);
      free(input);
      free(f0);
    }
    ms2_profile("analysis.load_or_generate_total", profile_started, filename);
    return chunk;
}
