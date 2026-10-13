/*
 * pi4decoder - a small, dependency-free PI4 beacon decoder.
 *
 * Copyright (C) 2026 Dawid Szymanski SQ6EMM
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PI4_H
#define PI4_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- protocol constants (OZ2M / DJ5HG "Next Generation Beacon" PI4) ---- */
#define PI4_NSYM      146        /* symbols per transmission                 */
#define PI4_NINFO     42         /* information bits (8 chars base 38)       */
#define PI4_NBITS     73         /* info + 31-bit zero tail                  */
#define PI4_NCHAR     8
#define PI4_SYMBOL_S  (1.0/6.0)  /* 166.667 ms; 360 symbols per minute       */
#define PI4_FS        12000      /* internal sample rate                     */
#define PI4_SYMLEN    2000       /* samples per symbol at PI4_FS             */

extern const char    pi4_alphabet[39];      /* "0123456789A..Z /"          */
extern const uint8_t pi4_sync[PI4_NSYM];    /* DJ5HG sync vector           */

typedef enum { PI4_V_PI4 = 0, PI4_V_80 = 1, PI4_V_96 = 2, PI4_V_120 = 3 } pi4_variant;
#define PI4_VARIANT_MASK_ALL 0x0f

/* ---- source coding / FEC (spec module) ---- */
int      pi4_pack(const char *msg, uint64_t *n);          /* 0 ok, -1 bad char/len */
void     pi4_unpack(uint64_t n, char out[PI4_NCHAR + 1]); /* trailing spaces trimmed */
void     pi4_encode_n(uint64_t n, uint8_t tones[PI4_NSYM]);
int      pi4_encode(const char *msg, uint8_t tones[PI4_NSYM]);
void     pi4_coded_bits(uint64_t n, uint8_t coded[PI4_NSYM]); /* codeword order */
extern const uint8_t pi4_interleave[PI4_NSYM];  /* codeword index -> channel symbol */
double   pi4_spacing_hz(pi4_variant v);          /* K * 12000 / 2048         */
int      pi4_variant_k(pi4_variant v);
const char *pi4_variant_name(pi4_variant v);
int      pi4_variant_parse(const char *s);       /* -1 if unknown            */

/* ---- sequential (Fano) decoder ----
 * metric[i][b]: integer branch metric for coded bit i (codeword order) taking
 * value b. Returns 1 if the search converged (bits[0..72] filled, info bits
 * first), 0 if the cycle budget ran out. *metric_out gets the final path
 * metric. */
int pi4_fano(const int16_t *metric, int delta, long maxcycles, int reverse,
             uint8_t bits[PI4_NBITS], long *metric_out, long *cycles_out);

/* ---- decoder ---- */
typedef struct {
    double center_hz;    /* CW-carrier convention: tone0 = center - spacing/2 (800) */
    double radius_hz;    /* tone0 searched +-radius around the convention (150) */
    double fmin_hz, fmax_hz; /* if fmax>fmin: explicit tone0 band instead (scan) */
    double t_lo_s, t_hi_s;   /* message start searched in [t_lo,t_hi] (file time) */
    unsigned variants;   /* bitmask of pi4_variant                            */
    int   max_spread;    /* widest spread window index (0..PI4_NSPREAD-1)     */
    int   max_cand;      /* candidate cells tried with the sequential decoder */
    long  fano_cycles;   /* cycle budget per attempt                          */
    double min_sync_z;   /* candidate gate on the sync statistic              */
    double min_fit_z;    /* decode gate (blind mode)                          */
    int   relax;         /* 1: accept also implausible texts at a higher gate */
    const char *const *calls; int ncalls; /* known-message list (0 = blind)   */
    double list_min_z;   /* gate for list mode                                */
    int   verbose;       /* debug prints to stderr                            */
    int   metric;        /* -1: automatic sequence; 0 rayleigh, 1 bessel, 2 gauss (testing) */
} pi4_params;

#define PI4_NSPREAD 6
extern const int pi4_spread_bins[PI4_NSPREAD];  /* 1,3,5,9,17,33 bins of fs/4096 */

typedef struct {
    char   text[PI4_NCHAR + 1];
    pi4_variant variant;
    double tone0_hz;      /* audio frequency of tone 0                         */
    double t_s;           /* message start, seconds from file start            */
    double snr_db;        /* S/N in 2500 Hz                                     */
    double sync_z;        /* sync-correlation z-score at the decoded cell      */
    double fit;           /* soft agreement of re-encoded message, -1..1       */
    double fit_z;         /* the same as a z-score                             */
    int    hard_errors;   /* symbols where the re-encoded data bit disagrees   */
    double spread_hz;     /* energy window used per tone                       */
    int    plausible;     /* text looks like a callsign / status message       */
    int    listed;        /* 1: found by the known-message matcher             */
    long   fano_metric;
    double energy;        /* signal energy explained by the decode (noise units) */
} pi4_result;

void pi4_params_default(pi4_params *p);

/* Whitened spectrogram: the decoder's working data. Built from audio, can be
 * accumulated over several minutes (averaging) before decoding. */
typedef struct pi4_spec pi4_spec;
pi4_spec *pi4_spec_build(const float *audio, size_t n, int fs, const pi4_params *p);
int       pi4_spec_add(pi4_spec *acc, const pi4_spec *s);   /* same geometry */
void      pi4_spec_free(pi4_spec *s);

/* Decode: audio (any of 6/8/12/24/48 kHz, mono float) -> results.
 * Returns the number of results (<= max_out). */
int pi4_decode(const float *audio, size_t n, int fs, const pi4_params *p,
               pi4_result *out, int max_out);
/* Decode from a (possibly averaged) spectrogram. audio may be NULL (then no
 * time-domain refinement is done; fs/n describe the audio if given). */
int pi4_decode_spec(const pi4_spec *s, const float *audio, size_t n, int fs,
                    const pi4_params *p, pi4_result *out, int max_out);

/* Resample mono audio to 12 kHz. Returns malloc'd buffer, *n_out samples. */
float *pi4_resample12k(const float *in, size_t n, int fs, size_t *n_out);

/* Text plausibility: 1 = callsign-like (optionally /suffix) or /status. */
int pi4_text_plausible(const char *text);

#ifdef __cplusplus
}
#endif
#endif
