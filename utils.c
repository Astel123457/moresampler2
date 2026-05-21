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
    opt->chanfreq = malloc(sizeof(FP_TYPE) * chanfreq_len);
    if (!opt->chanfreq) return -1;
    fread(opt->chanfreq, sizeof(FP_TYPE), chanfreq_len, f);
    fread(&opt->lip_radius, sizeof(FP_TYPE), 1, f);
    fread(&opt->f0_refine, sizeof(FP_TYPE), 1, f);
    fread(&opt->hm_method, sizeof(int), 1, f);
    fread(&opt->rel_winsize, sizeof(FP_TYPE), 1, f);
    return 0;
}

int save_llsm(llsm_chunk* chunk, const char* filename, llsm_aoptions* conf, int* fs, int* nbit) {
  // Detect which layer the chunk is in by inspecting the first frame.
  int* nfrm_p = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  int nfrm_val = nfrm_p ? *nfrm_p : 0;
  int layer = (nfrm_val > 0 && llsm_frame_checklayer1(chunk->frames[0])) ? 1 : 0;

  // For Layer 1 we need NSPEC from the chunk conf (set by llsm_chunk_tolayer1).
  int nspec = 0;
  if (layer == 1) {
    int* nspec_p = llsm_container_get(chunk->conf, LLSM_CONF_NSPEC);
    nspec = nspec_p ? *nspec_p : 0;
  }

  FILE* f = fopen(filename, "wb");
  if (!f) return -1;

  // Version 2 header: adds layer and nspec before nfrm.
  fwrite("LLSM2", 1, 5, f);
  int version = 2;
  fwrite(&version, sizeof(int), 1, f);
  fwrite(&layer,   sizeof(int), 1, f);
  fwrite(&nspec,   sizeof(int), 1, f);
  fwrite(&nfrm_val, sizeof(int), 1, f);
  fwrite(fs,   sizeof(int), 1, f);
  fwrite(nbit, sizeof(int), 1, f);
  write_conf(f, conf);

  // Frame data
  for (int i = 0; i < nfrm_val; ++i) {
    llsm_container* frame = chunk->frames[i];
    FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
    fwrite(f0, sizeof(FP_TYPE), 1, f);

    if (layer == 0) {
      // Layer 0: harmonic model (nhar, ampl, phse)
      llsm_hmframe* hm = llsm_container_get(frame, LLSM_FRAME_HM);
      fwrite(&hm->nhar, sizeof(int), 1, f);
      fwrite(hm->ampl,  sizeof(FP_TYPE), hm->nhar, f);
      fwrite(hm->phse,  sizeof(FP_TYPE), hm->nhar, f);
    } else {
      // Layer 1: Rd scalar, VTMAGN fparray, VSPHSE fparray
      FP_TYPE* rd = llsm_container_get(frame, LLSM_FRAME_RD);
      FP_TYPE rd_val = rd ? *rd : (FP_TYPE)1.0;
      fwrite(&rd_val, sizeof(FP_TYPE), 1, f);

      FP_TYPE* vtmagn = llsm_container_get(frame, LLSM_FRAME_VTMAGN);
      int vtmagn_n = (vtmagn && *f0 > 0) ? llsm_fparray_length(vtmagn) : 0;
      fwrite(&vtmagn_n, sizeof(int), 1, f);
      if (vtmagn_n > 0) fwrite(vtmagn, sizeof(FP_TYPE), vtmagn_n, f);

      FP_TYPE* vsphse = llsm_container_get(frame, LLSM_FRAME_VSPHSE);
      int vsphse_n = (vsphse && *f0 > 0) ? llsm_fparray_length(vsphse) : 0;
      fwrite(&vsphse_n, sizeof(int), 1, f);
      if (vsphse_n > 0) fwrite(vsphse, sizeof(FP_TYPE), vsphse_n, f);
    }

    // NM is identical in both layers
    llsm_nmframe* nm = llsm_container_get(frame, LLSM_FRAME_NM);
    fwrite(&nm->npsd, sizeof(int), 1, f);
    fwrite(nm->psd,   sizeof(FP_TYPE), nm->npsd, f);
    fwrite(&nm->nchannel, sizeof(int), 1, f);
    for (int j = 0; j < nm->nchannel; ++j) {
      fwrite(&nm->edc[j], sizeof(FP_TYPE), 1, f);
      llsm_hmframe* eenv = nm->eenv[j];
      fwrite(&eenv->nhar, sizeof(int), 1, f);
      fwrite(eenv->ampl,  sizeof(FP_TYPE), eenv->nhar, f);
      fwrite(eenv->phse,  sizeof(FP_TYPE), eenv->nhar, f);
    }
  }

  fclose(f);
  return 0;
}

llsm_chunk* read_llsm(const char* filename, int* nfrm, int* fs, int* nbit) {
  FILE* f = fopen(filename, "rb");
  if (!f) return NULL;

  char header[5];
  fread(header, 1, 5, f);
  if (strncmp(header, "LLSM2", 5) != 0) { fclose(f); return NULL; }

  int version;
  fread(&version, sizeof(int), 1, f);
  if (version != 1 && version != 2) { fclose(f); return NULL; }

  // Version 2 adds layer and nspec before nfrm.
  int layer = 0, nspec = 0;
  if (version == 2) {
    fread(&layer, sizeof(int), 1, f);
    fread(&nspec, sizeof(int), 1, f);
  }

  fread(nfrm, sizeof(int), 1, f);
  fread(fs,   sizeof(int), 1, f);
  fread(nbit, sizeof(int), 1, f);

  llsm_aoptions* aopt = llsm_create_aoptions();
  if (read_conf(f, aopt) != 0) {
    llsm_delete_aoptions(aopt);
    fclose(f);
    return NULL;
  }
  llsm_container* conf = llsm_aoptions_toconf(aopt, 44100.0 / 2);
  llsm_delete_aoptions(aopt);
  llsm_container_attach(conf, LLSM_CONF_NFRM,
    llsm_create_int(*nfrm), llsm_delete_int, llsm_copy_int);
  if (layer == 1 && nspec > 0) {
    llsm_container_attach(conf, LLSM_CONF_NSPEC,
      llsm_create_int(nspec), llsm_delete_int, llsm_copy_int);
  }
  llsm_chunk* chunk = llsm_create_chunk(conf, *nfrm);
  llsm_delete_container(conf);

  for (int i = 0; i < *nfrm; ++i) {
    llsm_container* frame = llsm_create_frame(0, 0, 0, 0);

    // F0 (same for both layers)
    FP_TYPE* f0 = malloc(sizeof(FP_TYPE));
    fread(f0, sizeof(FP_TYPE), 1, f);
    llsm_container_attach(frame, LLSM_FRAME_F0, f0, free, llsm_copy_fp);

    if (layer == 0) {
      // Layer 0: harmonic model
      int nhar;
      fread(&nhar, sizeof(int), 1, f);
      llsm_hmframe* hm = llsm_create_hmframe(nhar);
      fread(hm->ampl, sizeof(FP_TYPE), nhar, f);
      fread(hm->phse, sizeof(FP_TYPE), nhar, f);
      llsm_container_attach(frame, LLSM_FRAME_HM, hm, llsm_delete_hmframe, llsm_copy_hmframe);
    } else {
      // Layer 1: Rd, VTMAGN, VSPHSE
      // llsm_create_frame pre-attaches HM(nhar=0); remove it so
      // llsm_frame_checklayer0 won't mistake this for a Layer 0 frame.
      llsm_container_remove(frame, LLSM_FRAME_HM);
      FP_TYPE* rd = malloc(sizeof(FP_TYPE));
      fread(rd, sizeof(FP_TYPE), 1, f);
      llsm_container_attach(frame, LLSM_FRAME_RD, rd, llsm_delete_fp, llsm_copy_fp);

      int vtmagn_n;
      fread(&vtmagn_n, sizeof(int), 1, f);
      if (vtmagn_n > 0) {
        FP_TYPE* vtmagn = llsm_create_fparray(vtmagn_n);
        fread(vtmagn, sizeof(FP_TYPE), vtmagn_n, f);
        llsm_container_attach(frame, LLSM_FRAME_VTMAGN, vtmagn, llsm_delete_fparray, llsm_copy_fparray);
      }

      int vsphse_n;
      fread(&vsphse_n, sizeof(int), 1, f);
      if (vsphse_n > 0) {
        FP_TYPE* vsphse = llsm_create_fparray(vsphse_n);
        fread(vsphse, sizeof(FP_TYPE), vsphse_n, f);
        llsm_container_attach(frame, LLSM_FRAME_VSPHSE, vsphse, llsm_delete_fparray, llsm_copy_fparray);
      }
    }

    // NM (same for both layers)
    llsm_nmframe* nm = malloc(sizeof(llsm_nmframe));
    fread(&nm->npsd, sizeof(int), 1, f);
    nm->psd = malloc(sizeof(FP_TYPE) * nm->npsd);
    fread(nm->psd, sizeof(FP_TYPE), nm->npsd, f);
    fread(&nm->nchannel, sizeof(int), 1, f);
    nm->edc  = malloc(sizeof(FP_TYPE)      * nm->nchannel);
    nm->eenv = malloc(sizeof(llsm_hmframe*) * nm->nchannel);
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
    chunk->frames[i] = frame;
  }

  fclose(f);
  return chunk;
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

llsm_chunk* get_chunk_from_file(const char* filename, int* nfrm, int* fs, int* nbit, int* nx, int* nhop, int nosave, llsm_aoptions* opt_a, llsm_soptions* opt_s) {

    llsm_chunk* chunk = NULL;
    FP_TYPE* f0 = NULL;
    FP_TYPE* input = NULL;

    // Build expected .llsm2 path from input WAV path
    char* llsm_path = build_llsm_path(filename);
    if (!llsm_path) {
      printf("WARNING: Failed to build LLSM2 filename path. Analysis will proceed without caching.\n");
    }
    FILE* llsm_file = fopen(llsm_path, "rb");

    if (llsm_file) {
      // File exists — use cached analysis
      fclose(llsm_file);
      printf("Loading cached LLSM analysis: %s\n", llsm_path);
      chunk = read_llsm(llsm_path, nfrm, fs, nbit);
      free(llsm_path);
      if (!chunk) {
        printf("Failed to read .llsm2 file\n");
        return NULL;
      }
    } else {
      // No cache — analyze audio
      printf("Reading input WAV: %s\n", filename);
      input = wavread(filename, fs, nbit, nx);
      if (!input) { free(llsm_path); return NULL; }
    
      printf("Estimating F0\n");
      pyin_config param = pyin_init(*nhop);
      param.fmin = 50.0f;
      param.fmax = 800.0f;
      param.trange = 24;
      param.bias = 2;
      param.nf = ceil(*fs * 0.025);
      f0 = pyin_analyze(param, input, *nx, *fs, nfrm);
      if (!f0) { free(input); free(llsm_path); return NULL; }
    
      opt_a->thop = (FP_TYPE)(*nhop) / (*fs);
      opt_a->f0_refine = 1;
      opt_a->hm_method = LLSM_AOPTION_HMCZT;
    
      printf("Analysis\n");
      chunk = llsm_analyze(opt_a, input, *nx, *fs, f0, *nfrm, NULL);
      if (!chunk) { free(input); free(f0); free(llsm_path); return NULL; }
      if(llsm_path != NULL) {
        printf("Saving analysis result to cache: %s\n", llsm_path);
          if (!nosave) {
            if (save_llsm(chunk, llsm_path, opt_a, fs, nbit) != 0) {
              printf("Failed to save .llsm2 file.\n");
            }
          }
      }
      free(llsm_path);
      free(input);
      free(f0);
    }
    return chunk;
}