#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <ctype.h>
#include <ciglet/ciglet.h>
#include <libllsm/llsm.h>
#include <libpyin/pyin.h>

const char* version = "0.2.0";

// circular interpolation of two radian values
static FP_TYPE linterpc(FP_TYPE a, FP_TYPE b, FP_TYPE ratio) {
  FP_TYPE ax = cos_2(a);
  FP_TYPE ay = sin_2(a);
  FP_TYPE bx = cos_2(b);
  FP_TYPE by = sin_2(b);
  FP_TYPE cx = linterp(ax, bx, ratio);
  FP_TYPE cy = linterp(ay, by, ratio);
  return atan2(cy, cx);
}

static void interp_nmframe(llsm_nmframe* dst, llsm_nmframe* src,
  FP_TYPE ratio, int dst_voiced, int src_voiced) {
  for(int i = 0; i < dst -> npsd; i ++)
    dst -> psd[i] = linterp(dst -> psd[i], src -> psd[i], ratio);

  for(int b = 0; b < dst -> nchannel; b ++) {
    llsm_hmframe* srceenv = src -> eenv[b];
    llsm_hmframe* dsteenv = dst -> eenv[b];
    dst -> edc[b] = linterp(dst -> edc[b], src -> edc[b], ratio);
    int b_minnhar = min(srceenv -> nhar, dsteenv -> nhar);
    int b_maxnhar = max(srceenv -> nhar, dsteenv -> nhar);
    if(dsteenv -> nhar < b_maxnhar) {
      dsteenv -> ampl = realloc(dsteenv -> ampl, sizeof(FP_TYPE) * b_maxnhar);
      dsteenv -> phse = realloc(dsteenv -> phse, sizeof(FP_TYPE) * b_maxnhar);
    }
    for(int i = 0; i < b_minnhar; i ++) {
      dsteenv -> ampl[i] =
        linterp(dsteenv -> ampl[i], srceenv -> ampl[i], ratio);
      dsteenv -> phse[i] =
        linterpc(dsteenv -> phse[i], srceenv -> phse[i], ratio);
    }
    if(b_maxnhar == srceenv -> nhar) {
      for(int i = b_minnhar; i < b_maxnhar; i ++) {
        dsteenv -> ampl[i] = srceenv -> ampl[i];
        dsteenv -> phse[i] = srceenv -> phse[i];
      }
    }
    dsteenv -> nhar = b_maxnhar;
  }
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
  FILE* f = fopen(filename, "wb");
  if (!f) return -1;

  // Header
  fwrite("LLSM2", 1, 5, f);
  int version = 1;
  fwrite(&version, sizeof(int), 1, f);

  // Frame count
  int* nfrm = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
  fwrite(nfrm, sizeof(int), 1, f);
  fwrite(fs, sizeof(int), 1, f);
  fwrite(nbit, sizeof(int), 1, f);

  // Write chunk->conf
  write_conf(f, conf);

  // Frame data
  for (int i = 0; i < *nfrm; ++i) {
    llsm_container* frame = chunk->frames[i];

    // f0
    FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
    fwrite(f0, sizeof(FP_TYPE), 1, f);

    // HM Frame
    llsm_hmframe* hm = llsm_container_get(frame, LLSM_FRAME_HM);
    fwrite(&hm->nhar, sizeof(int), 1, f);
    fwrite(hm->ampl, sizeof(FP_TYPE), hm->nhar, f);
    fwrite(hm->phse, sizeof(FP_TYPE), hm->nhar, f);

    // NM Frame
    llsm_nmframe* nm = llsm_container_get(frame, LLSM_FRAME_NM);
    fwrite(&nm->npsd, sizeof(int), 1, f);
    fwrite(nm->psd, sizeof(FP_TYPE), nm->npsd, f);

    fwrite(&nm->nchannel, sizeof(int), 1, f);
    for (int j = 0; j < nm->nchannel; ++j) {
      fwrite(&nm->edc[j], sizeof(FP_TYPE), 1, f);

      llsm_hmframe* eenv = nm->eenv[j];
      fwrite(&eenv->nhar, sizeof(int), 1, f);
      fwrite(eenv->ampl, sizeof(FP_TYPE), eenv->nhar, f);
      fwrite(eenv->phse, sizeof(FP_TYPE), eenv->nhar, f);
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
  if (version != 1) { fclose(f); return NULL; }
  fread(nfrm, sizeof(int), 1, f);
  fread(fs, sizeof(int), 1, f);
  fread(nbit, sizeof(int), 1, f);
  // Read conf
  llsm_aoptions* aopt = llsm_create_aoptions();
  int conf_r = read_conf(f, aopt);
  if (conf_r != 0) {
    llsm_delete_aoptions(aopt);
    fclose(f);
    return NULL;
  }
  llsm_container* conf = llsm_aoptions_toconf(aopt, 44100.0 / 2);
  llsm_container_attach(conf, LLSM_CONF_NFRM,
    llsm_create_int(*nfrm), llsm_delete_int, llsm_copy_int);
  llsm_chunk* chunk = llsm_create_chunk(conf, *nfrm);
  

  for (int i = 0; i < *nfrm; ++i) {
    llsm_container* frame = llsm_create_frame(0, 0, 0, 0); // dummy init

    // f0
    FP_TYPE* f0 = malloc(sizeof(FP_TYPE));
    fread(f0, sizeof(FP_TYPE), 1, f);
    llsm_container_attach(frame, LLSM_FRAME_F0, f0, free, llsm_copy_fp);

    // HM
    int nhar;
    fread(&nhar, sizeof(int), 1, f);
    llsm_hmframe* hm = llsm_create_hmframe(nhar);
    fread(hm->ampl, sizeof(FP_TYPE), nhar, f);
    fread(hm->phse, sizeof(FP_TYPE), nhar, f);
    llsm_container_attach(frame, LLSM_FRAME_HM, hm, llsm_delete_hmframe, llsm_copy_hmframe);

    // NM
    llsm_nmframe* nm = malloc(sizeof(llsm_nmframe));
    fread(&nm->npsd, sizeof(int), 1, f);
    nm->psd = malloc(sizeof(FP_TYPE) * nm->npsd);
    fread(nm->psd, sizeof(FP_TYPE), nm->npsd, f);

    fread(&nm->nchannel, sizeof(int), 1, f);
    nm->edc = malloc(sizeof(FP_TYPE) * nm->nchannel);
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

#define LOG2DB (20.0 / 2.3025851)
#define mag2db(x) (log_2(x) * LOG2DB)

// dst <- (dst &> src)
static void interp_llsm_frame(llsm_container* dst, llsm_container* src,
  FP_TYPE ratio) {
# define EPS 1e-8
  FP_TYPE dst_f0 = *((FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_F0));
  FP_TYPE src_f0 = *((FP_TYPE*)llsm_container_get(src, LLSM_FRAME_F0));
  llsm_nmframe* dst_nm = llsm_container_get(dst, LLSM_FRAME_NM);
  llsm_nmframe* src_nm = llsm_container_get(src, LLSM_FRAME_NM);
  FP_TYPE* src_rd = llsm_container_get(src, LLSM_FRAME_RD);
  FP_TYPE* dst_rd = llsm_container_get(dst, LLSM_FRAME_RD);
  FP_TYPE* dst_vsphse = llsm_container_get(dst, LLSM_FRAME_VSPHSE);
  FP_TYPE* src_vsphse = llsm_container_get(src, LLSM_FRAME_VSPHSE);
  FP_TYPE* dst_vtmagn = llsm_container_get(dst, LLSM_FRAME_VTMAGN);
  FP_TYPE* src_vtmagn = llsm_container_get(src, LLSM_FRAME_VTMAGN);

  // always take the frequency of the voiced frame
  llsm_container* voiced = dst_f0 <= 0 && src_f0 <= 0 ? NULL :
    (src_f0 > 0 ? src : dst);
  int bothvoiced = dst_f0 > 0 && src_f0 > 0;

  int dstnhar = dst_vsphse == NULL ? 0 : llsm_fparray_length(dst_vsphse);
  int srcnhar = src_vsphse == NULL ? 0 : llsm_fparray_length(src_vsphse);
  int maxnhar = max(dstnhar, srcnhar);
  int minnhar = min(dstnhar, srcnhar);

  if(! bothvoiced && voiced == src) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(src_f0),
      llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(*src_rd),
      llsm_delete_fp, llsm_copy_fp);
  } else
  if(voiced == NULL) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(0),
      llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(1.0),
      llsm_delete_fp, llsm_copy_fp);
  }
  int nspec = dst_vtmagn != NULL ? llsm_fparray_length(dst_vtmagn) :
    (src_vtmagn != NULL ? llsm_fparray_length(src_vtmagn) : 0);

  if(bothvoiced) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(
      linterp(dst_f0, src_f0, ratio)), llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(
      linterp(*dst_rd, *src_rd, ratio)), llsm_delete_fp, llsm_copy_fp);

    FP_TYPE* vsphse = llsm_create_fparray(maxnhar);
    FP_TYPE* vtmagn = llsm_create_fparray(nspec);
    for(int i = 0; i < minnhar; i ++)
      vsphse[i] = linterpc(dst_vsphse[i], src_vsphse[i], ratio);
    for(int i = 0; i < nspec; i ++)
      vtmagn[i] = linterp(dst_vtmagn[i], src_vtmagn[i], ratio);
    if(dstnhar < srcnhar)
      for(int i = minnhar; i < maxnhar; i ++)
        vsphse[i] = src_vsphse[i];

    dst_vsphse = vsphse;
    dst_vtmagn = vtmagn;
    llsm_container_attach(dst, LLSM_FRAME_VSPHSE, dst_vsphse,
      llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(dst, LLSM_FRAME_VTMAGN, dst_vtmagn,
      llsm_delete_fparray, llsm_copy_fparray);
  } else if(voiced == src) {
    dst_vsphse = llsm_copy_fparray(src_vsphse);
    dst_vtmagn = llsm_copy_fparray(src_vtmagn);
    llsm_container_attach(dst, LLSM_FRAME_VSPHSE, dst_vsphse,
      llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(dst, LLSM_FRAME_VTMAGN, dst_vtmagn,
      llsm_delete_fparray, llsm_copy_fparray);
    FP_TYPE fade = mag2db(max(EPS, ratio));
    for(int i = 0; i < nspec; i ++) dst_vtmagn[i] += fade;
  } else {
    FP_TYPE fade = mag2db(max(EPS, 1.0 - ratio));
    for(int i = 0; i < nspec; i ++) dst_vtmagn[i] += fade;
  }
  for(int i = 0; i < nspec; i ++) dst_vtmagn[i] = max(-80, dst_vtmagn[i]);

  interp_nmframe(dst_nm, src_nm, ratio, dst_f0 > 0, src_f0 > 0);
# undef EPS
}

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

int getF0Contour(char *input, double *output)
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

void convert_cents_to_hz_offset(double* cents, int cents_len,
                                int nfrm, int nhop, int fs,
                                float* out_hz) {
    const float frame_duration_sec = (float)nhop / fs;
    const float input_interval_sec = 0.005f;  // 5ms per input offset

    for (int i = 0; i < nfrm; ++i) {
        float time_sec = i * frame_duration_sec;
        float idx = time_sec / input_interval_sec;

        int i0 = (int)idx;
        int i1 = i0 + 1;
        if (i1 >= cents_len) i1 = cents_len - 1;

        float frac = idx - i0;
        float cents_interp = cents[i0] * (1.0f - frac) + cents[i1] * frac;

        // Convert cents offset into Hz offset from base (0 cents = 0 Hz offset)
        float ratio = pow(2.0f, cents_interp / 1200.0f);
        out_hz[i] = (ratio - 1.0f);  // This is the *offset* in Hz (if applied to a base freq)
    }
}


// As a note, while experimenting i thought that this was the volume of the whole frame, not the vocal tract
// so this basically edits the tension of the vocal tract, not the volume of the whole frame
void llsm_normalize_volume(llsm_chunk* chunk, FP_TYPE target_peak) {
    int* total_frames = llsm_container_get(chunk->conf, LLSM_CONF_NFRM);
    FP_TYPE global_peak = 0.0f;

    // Step 1: Find the global peak magnitude (linear)
    for (int i = 0; i < *total_frames; ++i) {
        FP_TYPE* vt_magn = llsm_container_get(chunk->frames[i], LLSM_FRAME_VTMAGN);
        if (!vt_magn) continue;

        int nspec = llsm_fparray_length(vt_magn);
        for (int j = 0; j < nspec; ++j) {
            if (vt_magn[j] > global_peak)
                global_peak = vt_magn[j];
        }
    }

    // Step 2: Avoid divide-by-zero
    if (global_peak < 1e-9f || target_peak <= 0.0f)
        return;

    // Step 3: Compute gain in dB
    FP_TYPE linear_gain = target_peak / global_peak;
    FP_TYPE gain_db = 20.0f * log10(linear_gain);

    // Step 4: Apply dB gain to all magnitudes
    for (int i = 0; i < *total_frames; ++i) {
        FP_TYPE* vt_magn = llsm_container_get(chunk->frames[i], LLSM_FRAME_VTMAGN);
        if (!vt_magn) continue;

        int nspec = llsm_fparray_length(vt_magn);
        for (int j = 0; j < nspec; ++j) {
            vt_magn[j] += gain_db;
        }
    }
}

void normalize_waveform(FP_TYPE* waveform, int length, FP_TYPE target_peak) {
    FP_TYPE peak = 0.0f;
    for (int i = 0; i < length; ++i) {
        FP_TYPE abs_val = fabsf(waveform[i]);
        if (abs_val > peak) peak = abs_val;
    }

    if (peak < 1e-9f) return; // avoid divide-by-zero

    FP_TYPE scale = target_peak / peak;
    for (int i = 0; i < length; ++i)
        waveform[i] *= scale;
}

int parse_tempo(const char *tempo_str) {
    if (tempo_str[0] == '!') {
        tempo_str++;  // skip the '!'
    }
    return atoi(tempo_str);
}

typedef struct {
    char* input;  // Path to the input audio file
    char* output; // Path to the output audio file
    float tone; // Musical pitch of the note to be resampled, in hertz
    int velocity; // velocity of the consontant area
    char* flags; // raw string of resampler flags
    float offset; // offset in milliseconds
    float length; // length of the note in milliseconds
    float consonant; // length of the consonant area in milliseconds
    float cutoff; // cutoff frequency in hertz
    int volume; // volume of the note, 0-100
    int modulation; // modulation value, 0-100
    int tempo; // tempo in beats per minute
    char* pitch_curve; // pitch curve data
} resampler_data;

int resample(resampler_data* data) {
  // Allocate and load pitch curve
  double* f0_curve = malloc(sizeof(double) * 3000);
  if (!f0_curve) return 1;
  int pit_len = getF0Contour(data->pitch_curve, f0_curve);
  if (!pit_len) { free(f0_curve); return 1; }

  // Build expected .llsm2 path from input WAV path
  char llsm_path[1024];
  snprintf(llsm_path, sizeof(llsm_path), "%s", data->input);
  char* ext = strrchr(llsm_path, '.');
  if (ext) strcpy(ext, ".llsm2");  // Replace extension

  // Check for existing .llsm2 (ignore .llsm)
  FILE* llsm_file = fopen(llsm_path, "rb");

  llsm_aoptions* opt_a = llsm_create_aoptions();
  llsm_soptions* opt_s = llsm_create_soptions(44100.0f);
  llsm_chunk* chunk = NULL;
  int nhop = 128;
  int fs = 0, nbit = 0, nx = 0;
  float* input = NULL;
  FP_TYPE* f0 = NULL;
  int nfrm = 0;

  if (llsm_file) {
    // File exists — use cached analysis
    fclose(llsm_file);
    printf("Loading cached LLSM analysis: %s\n", llsm_path);
    chunk = read_llsm(llsm_path, &nfrm, &fs, &nbit);
    
    if (!chunk) {
      printf("Failed to read .llsm2 file\n");
      free(f0_curve);
      return 1;
    }
  } else {
    // No cache — analyze audio
    printf("Reading input WAV: %s\n", data->input);
    input = wavread(data->input, &fs, &nbit, &nx);
    if (!input) { free(f0_curve); return 1; }

    printf("Estimating F0\n");
    pyin_config param = pyin_init(nhop);
    param.fmin = 50.0f;
    param.fmax = 800.0f;
    param.trange = 24;
    param.bias = 2;
    param.nf = ceil(fs * 0.025);
    f0 = pyin_analyze(param, input, nx, fs, &nfrm);
    if (!f0) { free(input); free(f0_curve); return 1; }

    opt_a->thop = (FP_TYPE)nhop / fs;
    opt_a->f0_refine = 1;
    opt_a->hm_method = LLSM_AOPTION_HMCZT;

    printf("Analysis\n");
    chunk = llsm_analyze(opt_a, input, nx, fs, f0, nfrm, NULL);
    if (!chunk) { free(input); free(f0); free(f0_curve); return 1; }

    printf("Saving analysis result to cache: %s\n", llsm_path);
    if (save_llsm(chunk, llsm_path, opt_a, &fs, &nbit) != 0) {
      printf("Failed to save .llsm2 file.\n");
    }

    free(input);
    free(f0);
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
  if (end_frame > nfrm) end_frame = nfrm;
  if (end_frame <= start_frame) end_frame = start_frame + 1;

  // Calculate consonant frames (unstretched)
  int consonant_frames = (int)round((data->consonant / 1000.0) * fs / nhop);
  if (consonant_frames > end_frame - start_frame)
      consonant_frames = end_frame - start_frame;
  int sample_frames = end_frame - start_frame;
  // Calculate total output frames to match data->length (ms)
  int total_frames = (int)round((data->length / 1000.0) * fs / nhop);
  if (total_frames < consonant_frames) total_frames = consonant_frames + 1;
  float f0_array[total_frames];
  convert_cents_to_hz_offset(f0_curve, pit_len, total_frames, nhop, fs, f0_array);
  llsm_container* conf_new = llsm_copy_container(chunk -> conf);
  llsm_container_attach(conf_new, LLSM_CONF_NFRM,
    llsm_create_int(total_frames), llsm_delete_int, llsm_copy_int);
  llsm_chunk* chunk_new = llsm_create_chunk(conf_new, total_frames);
  llsm_delete_container(conf_new);
  int no_stretch = 0;
  // Copy consonant area directly
  if (total_frames <= sample_frames) {
    for (int i = 0; i < total_frames; i++) {
      chunk_new->frames[i] = llsm_copy_container(chunk->frames[start_frame + i]);
    }
    no_stretch = 1;
  } else {
    for (int i = 0; i < sample_frames; i++) {
      chunk_new->frames[i] = llsm_copy_container(chunk->frames[start_frame + i]);
    }
  }
  llsm_chunk_tolayer1(chunk_new, 2048);
  llsm_chunk_phasepropagate(chunk_new, -1);
  // Loop the vowel area instead of stretching
  if (no_stretch == 0) {
    // Only stretch the vowel area (after consonant_frames)
    int vowel_sample_frames = sample_frames - consonant_frames;
    int vowel_total_frames = total_frames - consonant_frames;
    for (int i = consonant_frames; i < total_frames; i++) {
      // Map output frame i to input frame in the vowel area
      FP_TYPE mapped = (FP_TYPE)(i - consonant_frames) * vowel_sample_frames / vowel_total_frames;
      int base = consonant_frames + (int)mapped;
      FP_TYPE ratio = mapped - (int)mapped;
      int residx = base + rand() % 5 - 2;
      residx = max(consonant_frames, min(consonant_frames + vowel_sample_frames - 1, residx));
      base = min(base, consonant_frames + vowel_sample_frames - 2);
      chunk_new->frames[i] = llsm_copy_container(chunk_new->frames[base]);
      interp_llsm_frame(
        chunk_new->frames[i], chunk_new->frames[base + 1], ratio);
      FP_TYPE* resvec = llsm_container_get(chunk_new->frames[residx],
        LLSM_FRAME_PSDRES);
      if (resvec != NULL) {
        llsm_container_attach(chunk_new->frames[i], LLSM_FRAME_PSDRES,
          llsm_copy_fparray(resvec), llsm_delete_fparray, llsm_copy_fparray);
      }
    }
  }
  
  // Set all f0 to target tone
  for(int i = 0; i < total_frames; i ++) {
    llsm_container_attach(chunk_new->frames[i], LLSM_FRAME_HM, NULL, NULL, NULL);
    FP_TYPE* f0_i = llsm_container_get(chunk_new->frames[i], LLSM_FRAME_F0);
    FP_TYPE old_f0 = f0_i[0];
    if (f0_i[0] == 0.00f) {
      continue;
    }
    else {
      f0_i[0] = data->tone * (1.0 + f0_array[i]);
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
  llsm_normalize_volume(chunk_new, 20.00f);
  llsm_chunk_phasepropagate(chunk_new, 1);
  llsm_chunk_tolayer0(chunk_new);
  printf("Synthesis\n");
  
  llsm_output* out = llsm_synthesize(opt_s, chunk_new);

  if (!out || !out->y) {
      printf("Failed to synthesize output: %d\n", out);
      return 1;
  }

  normalize_waveform(out->y, out->ny, 0.60f);

  float scale = data->velocity / 100.0f;
  for (int i = 0; i < out->ny; ++i)
    out->y[i] *= scale;

  wavwrite(out->y, out->ny, fs, nbit, data->output);

  llsm_delete_output(out);
  llsm_delete_chunk(chunk);
  llsm_delete_chunk(chunk_new);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0_curve);
  return 0;
}

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
    if (argc == 14) { // user wants the resampler mode
        resampler_data data;
        data.input = argv[1];
        data.output = argv[2];
        data.tone = note_to_frequency(argv[3]); // note is passed as a string (A4), so we need to convert it
        data.velocity = atoi(argv[4]);
        data.flags = argv[5]; // flags are passed as a string, e.g. "Mt50"
        data.offset = atof(argv[6]);
        data.length = atof(argv[7]);
        data.consonant = atof(argv[8]);
        data.cutoff = atof(argv[9]);
        data.volume = atoi(argv[10]);
        data.modulation = atoi(argv[11]);
        data.tempo = parse_tempo(argv[12]); //since tempo has a special format, we need to parse it
        data.pitch_curve = argv[13]; // pitch curve data as a string
        return resample(&data);
    }
    printf("Invalid arguments. Expected 14 arguments, got %d.\n", argc);
    return 0;
}