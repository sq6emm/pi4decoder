/*
 * PI4 decoder core: whitened spectrogram -> sync search over
 * (time, tone-0 frequency, tone spacing variant, spread window) ->
 * refinement -> soft metrics -> sequential decoding -> gates.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "pi4.h"
#include "pi4_fft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define HOP      250                      /* 1/8 symbol at 12 kHz            */
#define SPH      8                        /* hops per symbol                 */
#define NFFT     4096                     /* 2000 samples zero padded        */
#define BIN_HZ   (12000.0 / NFFT)         /* 2.9297 Hz                       */
#define MAX_CAND 64
#define CAND_PER_HYP 12

const int pi4_spread_bins[PI4_NSPREAD] = { 1, 3, 5, 9, 17, 33 };

struct pi4_spec {
    int nhop, h0;        /* columns; column c is file hop h0+c               */
    int b0, nb;          /* first bin, bin count                             */
    int navg;
    float *p;            /* whitened power [nhop][nb]                        */
    float *noise;        /* raw noise power per bin [nb]                     */
    int hvalid0, hvalid1;/* columns fully inside the audio                   */
};

void pi4_params_default(pi4_params *p)
{
    memset(p, 0, sizeof *p);
    p->center_hz = 800.0;
    p->radius_hz = 150.0;
    p->t_lo_s = -2.5; p->t_hi_s = 2.5;
    p->variants = PI4_VARIANT_MASK_ALL;
    p->max_spread = PI4_NSPREAD - 1;
    p->max_cand = 12;
    p->fano_cycles = 200000;
    p->min_sync_z = 4.0;
    p->min_fit_z = 6.0;
    p->list_min_z = 6.5;
    p->metric = -1;
}

/* ------------------------------------------------------------------ */
/* resampling to 12 kHz: windowed-sinc polyphase, rational L/M          */

static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

float *pi4_resample12k(const float *in, size_t n, int fs, size_t *n_out)
{
    const int H = 24;                 /* half kernel width in input samples */
    int g = gcd(fs, PI4_FS), L = PI4_FS / g, M = fs / g, ph, k;
    double c = fs > PI4_FS ? (double)PI4_FS / fs : 1.0;
    size_t nout = (size_t)((double)n * L / M), m;
    float *out = (float *)malloc(sizeof(float) * (nout + 1));
    float *ker;
    if (!out) return NULL;
    if (L == 1 && M == 1) { memcpy(out, in, sizeof(float) * n); *n_out = n; return out; }
    c *= 0.92;
    ker = (float *)malloc(sizeof(float) * (size_t)L * (2 * H + 1));
    for (ph = 0; ph < L; ph++) {
        double frac = (double)ph / L;             /* t = i0 + frac */
        for (k = -H; k <= H; k++) {
            double u = k - frac, w, s;            /* j - t with j = i0 + k */
            if (fabs(u) < 1e-9) s = c;
            else s = c * sin(M_PI * c * u) / (M_PI * c * u);
            w = 0.5 * (1 + cos(M_PI * u / (H + 1)));
            ker[ph * (2 * H + 1) + k + H] = (float)(s * w);
        }
    }
    for (m = 0; m < nout; m++) {
        size_t num = m * (size_t)M;
        long i0 = (long)(num / (size_t)L);
        const float *kp = ker + (num % (size_t)L) * (2 * H + 1);
        double acc = 0;
        for (k = -H; k <= H; k++) {
            long j = i0 + k;
            if (j >= 0 && (size_t)j < n) acc += in[j] * kp[k + H];
        }
        out[m] = (float)acc;
    }
    free(ker);
    *n_out = nout;
    return out;
}

/* ------------------------------------------------------------------ */
/* spectrogram                                                          */

static void band_for(const pi4_params *p, double *fmin, double *fmax, int *bt0, int *bt1)
{
    /* tone0 search band [bt0,bt1) in bins and the full tone band [fmin,fmax] */
    double lo = 1e9, hi = -1e9, tlo = 1e9, thi = -1e9;
    int v;
    for (v = 0; v < 4; v++) {
        double sp, t0lo, t0hi;
        if (!(p->variants & (1u << v))) continue;
        sp = pi4_spacing_hz((pi4_variant)v);
        if (p->fmax_hz > p->fmin_hz) { t0lo = p->fmin_hz; t0hi = p->fmax_hz; }
        else { t0lo = p->center_hz - sp / 2 - p->radius_hz; t0hi = p->center_hz - sp / 2 + p->radius_hz; }
        if (t0lo < tlo) tlo = t0lo;
        if (t0hi > thi) thi = t0hi;
        if (t0lo < lo) lo = t0lo;
        if (t0hi + 3 * sp > hi) hi = t0hi + 3 * sp;
    }
    *fmin = lo; *fmax = hi;
    *bt0 = (int)floor(tlo / BIN_HZ); *bt1 = (int)ceil(thi / BIN_HZ) + 1;
}

/* 1 when tone 0 lies outside the requested search band of this variant (by
 * more than half a bin): a signal caught at the band edge and re-centred.
 * It is kept and flagged, not dropped. */
static int outside_band(const pi4_params *p, pi4_variant v, double f0)
{
    double sp = pi4_spacing_hz(v), lo, hi;
    if (p->fmax_hz > p->fmin_hz) { lo = p->fmin_hz; hi = p->fmax_hz; }
    else { lo = p->center_hz - sp / 2 - p->radius_hz; hi = p->center_hz - sp / 2 + p->radius_hz; }
    return f0 < lo - BIN_HZ / 2 || f0 > hi + BIN_HZ / 2;
}

static int cmp_float(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return x < y ? -1 : x > y;
}

pi4_spec *pi4_spec_build(const float *audio, size_t n, int fs, const pi4_params *p)
{
    pi4_spec *s;
    pi4_fft fft;
    double fmin, fmax;
    int bt0, bt1, margin = pi4_spread_bins[PI4_NSPREAD - 1] / 2 + 2;
    int hlo, hhi, h, b;
    float *buf, *out, *col;
    const float *x = audio;
    float *own = NULL;
    size_t nx = n;

    if (fs != PI4_FS) {
        own = pi4_resample12k(audio, n, fs, &nx);
        if (!own) return NULL;
        x = own;
    }
    s = (pi4_spec *)calloc(1, sizeof *s);
    band_for(p, &fmin, &fmax, &bt0, &bt1);
    s->b0 = (int)floor(fmin / BIN_HZ) - margin;
    if (s->b0 < 1) s->b0 = 1;
    s->nb = (int)ceil(fmax / BIN_HZ) + margin + 1 - s->b0;
    if (s->b0 + s->nb > NFFT / 2) s->nb = NFFT / 2 - s->b0;
    hlo = (int)floor(p->t_lo_s * PI4_FS / HOP) - 2;
    hhi = (int)ceil(p->t_hi_s * PI4_FS / HOP) + 2 + PI4_NSYM * SPH;
    if (hlo < 0) hlo = 0;
    if (hhi < hlo + PI4_NSYM * SPH + 1) hhi = hlo + PI4_NSYM * SPH + 1;
    s->h0 = hlo; s->nhop = hhi - hlo + 1;
    s->navg = 1;
    s->p = (float *)calloc((size_t)s->nhop * (size_t)s->nb, sizeof(float));
    s->noise = (float *)calloc((size_t)s->nb, sizeof(float));
    if (pi4_fft_init(&fft, NFFT / 2)) { free(own); return NULL; }
    buf = (float *)malloc(sizeof(float) * NFFT);
    out = (float *)malloc(sizeof(float) * (NFFT + 2));
    s->hvalid0 = s->nhop; s->hvalid1 = 0;
    for (h = 0; h < s->nhop; h++) {
        long start = (long)(s->h0 + h) * HOP;
        int i;
        if (start >= (long)nx) continue;
        if (start + PI4_SYMLEN <= (long)nx) {
            if (h < s->hvalid0) s->hvalid0 = h;
            s->hvalid1 = h + 1;
        }
        memset(buf, 0, sizeof(float) * NFFT);
        for (i = 0; i < PI4_SYMLEN && start + i < (long)nx; i++) buf[i] = x[start + i];
        pi4_fft_real(&fft, buf, out);
        col = s->p + (size_t)h * s->nb;
        for (b = 0; b < s->nb; b++) {
            float re = out[2 * (s->b0 + b)], im = out[2 * (s->b0 + b) + 1];
            col[b] = re * re + im * im;
        }
    }
    free(buf); free(out); pi4_fft_free(&fft);
    free(own);

    /* noise floor: per-bin median over valid columns, smoothed over bins */
    {
        int nv = s->hvalid1 - s->hvalid0, k;
        float *tmp = (float *)malloc(sizeof(float) * (size_t)(nv > 0 ? nv : 1));
        float *med = (float *)malloc(sizeof(float) * (size_t)s->nb);
        for (b = 0; b < s->nb; b++) {
            if (nv <= 0) { med[b] = 1; continue; }
            for (h = 0; h < nv; h++) tmp[h] = s->p[(size_t)(s->hvalid0 + h) * s->nb + b];
            qsort(tmp, (size_t)nv, sizeof(float), cmp_float);
            /* 30th percentile of an exponential is 0.3567 x mean; robust to a
             * strong signal occupying up to ~70% of the hops in a tone bin */
            med[b] = tmp[(int)(nv * 0.30)] / 0.356675f;
        }
        for (b = 0; b < s->nb; b++) {             /* running median over 13 bins: a strong
                                                     tone's main lobe biases <= 5 of them */
            float win[13]; int cnt = 0;
            for (k = b - 6; k <= b + 6; k++) if (k >= 0 && k < s->nb) win[cnt++] = med[k];
            qsort(win, (size_t)cnt, sizeof(float), cmp_float);
            s->noise[b] = win[cnt / 2];
            if (s->noise[b] <= 0) s->noise[b] = 1e-20f;
        }
        for (h = 0; h < s->nhop; h++) {
            col = s->p + (size_t)h * s->nb;
            for (b = 0; b < s->nb; b++) {
                col[b] /= s->noise[b];
                if (col[b] > 1e4f) col[b] = 1e4f;      /* cap at 40 dB over the floor */
            }
        }
        if (p->verbose) {
            double m = 0, v = 0; size_t cnt = 0;
            for (h = s->hvalid0; h < s->hvalid1; h++)
                for (b = 0; b < s->nb; b++) { float q = s->p[(size_t)h * s->nb + b]; m += q; v += (double)q * q; cnt++; }
            m /= (double)cnt; v = v / (double)cnt - m * m;
            fprintf(stderr, "spectrogram: %d hops x %d bins (bin0 %d), valid %d..%d, whitened mean %.3f var %.3f\n",
                    s->nhop, s->nb, s->b0, s->hvalid0, s->hvalid1, m, v);
        }
        free(tmp); free(med);
    }
    return s;
}

int pi4_spec_add(pi4_spec *acc, const pi4_spec *s)
{
    size_t i, n;
    if (acc->nhop != s->nhop || acc->nb != s->nb || acc->b0 != s->b0 || acc->h0 != s->h0) return -1;
    n = (size_t)acc->nhop * acc->nb;
    for (i = 0; i < n; i++) acc->p[i] += s->p[i];
    if (s->hvalid0 < acc->hvalid0) acc->hvalid0 = s->hvalid0;
    if (s->hvalid1 > acc->hvalid1) acc->hvalid1 = s->hvalid1;
    acc->navg += s->navg;
    return 0;
}

void pi4_spec_free(pi4_spec *s)
{
    if (!s) return;
    free(s->p); free(s->noise); free(s);
}

/* ------------------------------------------------------------------ */
/* candidate cells                                                      */

typedef struct {
    int v, w;            /* variant, spread index                            */
    int h, b;            /* column, bin of tone 0 (spectrogram coordinates)  */
    double dh, db;       /* sub-cell refinement                              */
    double z;            /* sync z-score                                     */
    int call;            /* list mode: index of the matched message, -1      */
} cand;

/* robust z-normalisation statistics of a score array */
/* median and 1.4826*MAD of a (subsampled) score array: immune to a strong
 * signal or interferer occupying a sizeable fraction of the cells */
static void robust_stats(const float *x, size_t n, double *mean, double *sd)
{
    const size_t NMAX = 32768;
    size_t step = n > NMAX ? n / NMAX : 1, m = n / step, i;
    float *t = (float *)malloc(sizeof(float) * (m ? m : 1));
    double med;
    if (!t || m < 10) { *mean = 0; *sd = 1; free(t); return; }
    for (i = 0; i < m; i++) t[i] = x[i * step];
    qsort(t, m, sizeof(float), cmp_float);
    med = t[m / 2];
    for (i = 0; i < m; i++) t[i] = (float)fabs(t[i] - med);
    qsort(t, m, sizeof(float), cmp_float);
    *mean = med;
    *sd = 1.4826 * t[m / 2];
    if (*sd <= 0) *sd = 1;
    free(t);
}

static int cand_cmp(const void *a, const void *b)
{
    double x = ((const cand *)a)->z, y = ((const cand *)b)->z;
    return x > y ? -1 : x < y;
}

/* add a local maximum to the list (keeps the best MAX_CAND) */
static void push_cand(cand *list, int *n, const cand *c)
{
    int i;
    for (i = 0; i < *n; i++) {
        /* same cell (any variant), or the same frequency within 30 s under the
         * same hypothesis (a beacon sends once a minute; repeated peaks at one
         * frequency are an interferer): keep the best */
        if (list[i].call == c->call && abs(list[i].b - c->b) <= 4 &&
            (abs(list[i].h - c->h) <= SPH ||
             (list[i].v == c->v && list[i].w == c->w && abs(list[i].h - c->h) < 30 * PI4_FS / HOP))) {
            if (c->z > list[i].z) list[i] = *c;
            return;
        }
    }
    if (*n < MAX_CAND) { list[(*n)++] = *c; return; }
    {   int worst = 0;
        for (i = 1; i < *n; i++) if (list[i].z < list[worst].z) worst = i;
        if (c->z > list[worst].z) list[worst] = *c;
    }
}

/* Sync (call < 0) or known-message (call >= 0, tones given) correlation
 * search over one (variant, spread) hypothesis; local maxima above zmin go
 * to the candidate list. The spread window is applied row by row (a box sum
 * over w bins), so no copy of the spectrogram is needed. */
static void search_cells(const pi4_spec *s, const pi4_params *p,
                         int v, int w, int call, const uint8_t *tones, double zmin,
                         cand *list, int *ncand)
{
    int spb = 2 * pi4_variant_k((pi4_variant)v);          /* spacing in bins */
    int wb = pi4_spread_bins[w], r = wb / 2;
    int bt0, bt1, hc0, hc1, nbt, nhc, h, b, n, hmax;
    double fmin, fmax, mean, sd;
    float *row, *line, *score;
    cand *local = (cand *)malloc(sizeof(cand) * MAX_CAND); int nlocal = 0;
    band_for(p, &fmin, &fmax, &bt0, &bt1);
    bt0 -= s->b0; bt1 -= s->b0;
    if (bt0 < r + 1) bt0 = r + 1;                          /* no edge-clipped windows */
    if (bt1 > s->nb - 3 * spb - r) bt1 = s->nb - 3 * spb - r;
    nbt = bt1 - bt0;
    hc0 = (int)floor(p->t_lo_s * PI4_FS / HOP) - s->h0;
    hc1 = (int)ceil(p->t_hi_s * PI4_FS / HOP) - s->h0 + 1;
    if (hc0 < 0) hc0 = 0;
    if (hc1 > s->nhop - PI4_NSYM * SPH) hc1 = s->nhop - PI4_NSYM * SPH;
    nhc = hc1 - hc0;
    if (nbt <= 0 || nhc <= 0) { free(local); return; }

    row = (float *)malloc(sizeof(float) * (size_t)s->nb);
    line = (float *)malloc(sizeof(float) * (size_t)nbt);
    score = (float *)calloc((size_t)nhc * (size_t)nbt, sizeof(float));
    hmax = hc1 - 1 + SPH * (PI4_NSYM - 1);
    if (hmax >= s->nhop) hmax = s->nhop - 1;
    for (h = hc0; h <= hmax; h++) {
        const float *col = s->p + (size_t)h * s->nb;
        int nlo = (h - hc1 + 1 + SPH - 1) / SPH, nhi = (h - hc0) / SPH;
        if (nlo < 0) nlo = 0;
        if (nhi > PI4_NSYM - 1) nhi = PI4_NSYM - 1;
        if (nlo > nhi) continue;
        if (wb == 1) memcpy(row, col, sizeof(float) * (size_t)s->nb);
        else {                                             /* box sum over wb bins */
            double acc = 0;
            for (b = 0; b < r && b < s->nb; b++) acc += col[b];
            for (b = 0; b < s->nb; b++) {
                if (b + r < s->nb) acc += col[b + r];
                if (b - r - 1 >= 0) acc -= col[b - r - 1];
                row[b] = (float)acc;
            }
        }
        {
            const float *c0 = row + bt0;
            if (call < 0) {
                /* per-symbol discriminant: (even pair - odd pair) */
                for (b = 0; b < nbt; b++) line[b] = c0[b] + c0[b + 2 * spb] - c0[b + spb] - c0[b + 3 * spb];
                for (n = nlo; n <= nhi; n++) {
                    float *sc = score + (size_t)(h - SPH * n - hc0) * nbt;
                    if (pi4_sync[n]) for (b = 0; b < nbt; b++) sc[b] -= line[b];
                    else             for (b = 0; b < nbt; b++) sc[b] += line[b];
                }
            } else {
                /* matched filter: 4/3 P[tone] - 1/3 sum(all four) */
                for (b = 0; b < nbt; b++)
                    line[b] = -(1.0f / 3) * (c0[b] + c0[b + spb] + c0[b + 2 * spb] + c0[b + 3 * spb]);
                for (n = nlo; n <= nhi; n++) {
                    float *sc = score + (size_t)(h - SPH * n - hc0) * nbt;
                    const float *ct = c0 + tones[n] * spb;
                    for (b = 0; b < nbt; b++) sc[b] += line[b] + (4.0f / 3) * ct[b];
                }
            }
        }
    }
    free(row); free(line);
    robust_stats(score, (size_t)nhc * (size_t)nbt, &mean, &sd);
    for (h = 0; h < nhc; h++) {
        const float *sc = score + (size_t)h * nbt;
        for (b = 0; b < nbt; b++) {
            double z = (sc[b] - mean) / sd, zl, zr, zu, zd;
            cand c;
            if (z < zmin) continue;
            zl = b > 0 ? sc[b - 1] : -1e30; zr = b + 1 < nbt ? sc[b + 1] : -1e30;
            zu = h > 0 ? sc[b - nbt] : -1e30; zd = h + 1 < nhc ? sc[b + nbt] : -1e30;
            if (sc[b] < zl || sc[b] < zr || sc[b] < zu || sc[b] < zd) continue;
            c.v = v; c.w = w; c.h = hc0 + h; c.b = bt0 + b; c.z = z; c.call = call;
            /* parabolic sub-cell refinement */
            c.db = 0; c.dh = 0;
            if (zl > -1e29 && zr > -1e29) {
                double den = zl - 2 * sc[b] + zr;
                if (den < 0) c.db = 0.5 * (zl - zr) / den;
            }
            if (zu > -1e29 && zd > -1e29) {
                double den = zu - 2 * sc[b] + zd;
                if (den < 0) c.dh = 0.5 * (zu - zd) / den;
            }
            if (c.db > 0.5) c.db = 0.5;
            if (c.db < -0.5) c.db = -0.5;
            if (c.dh > 0.5) c.dh = 0.5;
            if (c.dh < -0.5) c.dh = -0.5;
            push_cand(local, &nlocal, &c);
        }
    }
    /* at most CAND_PER_HYP cells per hypothesis reach the global list, so one
     * strong interferer cannot crowd out a real signal under another one */
    qsort(local, (size_t)nlocal, sizeof(cand), cand_cmp);
    for (h = 0; h < nlocal && h < CAND_PER_HYP; h++) push_cand(list, ncand, &local[h]);
    free(score); free(local);
}

/* ------------------------------------------------------------------ */
/* per-symbol tone powers at a candidate                                */

/* Goertzel-free direct DFT of one symbol at 4 exact tone frequencies */
static void dft_tones(const float *x, size_t nx, long start, double f0, double sp,
                      float pw[4])
{
    int t, i;
    for (t = 0; t < 4; t++) {
        double f = f0 + t * sp, w = 2 * M_PI * f / PI4_FS;
        double cw = 2 * cos(w), s1 = 0, s2 = 0, s0;
        double re, im;
        for (i = 0; i < PI4_SYMLEN; i++) {
            long j = start + i;
            double v = (j >= 0 && (size_t)j < nx) ? x[j] : 0.0;
            s0 = v + cw * s1 - s2; s2 = s1; s1 = s0;
        }
        re = s1 - s2 * cos(w); im = s2 * sin(w);
        pw[t] = (float)(re * re + im * im);
    }
}

/* exact sync score for the fine search */
static double exact_sync(const pi4_spec *s, const float *x, size_t nx, long t0,
                         double f0, double sp, float (*pw)[4])
{
    double acc = 0; int n, t;
    for (n = 0; n < PI4_NSYM; n++) {
        float q[4];
        dft_tones(x, nx, t0 + (long)n * PI4_SYMLEN, f0, sp, q);
        for (t = 0; t < 4; t++) {
            int b = (int)floor((f0 + t * sp) / BIN_HZ + 0.5) - s->b0;
            if (b < 0) b = 0;
            if (b >= s->nb) b = s->nb - 1;
            q[t] /= s->noise[b];
            if (pw) pw[n][t] = q[t];
        }
        acc += pi4_sync[n] ? (q[1] + q[3] - q[0] - q[2]) : (q[0] + q[2] - q[1] - q[3]);
    }
    return acc;
}

/* self-normalised sync statistic of a set of per-symbol tone powers: the
 * sync-consistent pair minus the other pair, over the noise sd measured on
 * the pair that cannot carry signal */
static double sync_stat(float (*pw)[4])
{
    double acc = 0, var = 0; int n;
    for (n = 0; n < PI4_NSYM; n++) {
        int sy = pi4_sync[n];
        double e = pw[n][0] + pw[n][2], o = pw[n][1] + pw[n][3];
        double wd = pw[n][1 - sy] - pw[n][3 - sy];
        acc += sy ? o - e : e - o;
        var += wd * wd;
    }
    return var > 0 ? acc / sqrt(var) : 0;
}

/* tone powers from the spectrogram at an integer cell, summed over a
 * window of w bins around each tone */
static void spec_powers_w(const pi4_spec *s, int h, int b, int spb, int w, float (*pw)[4])
{
    int n, t, j, r = w / 2;
    for (n = 0; n < PI4_NSYM; n++) {
        const float *col = s->p + (size_t)(h + SPH * n) * s->nb;
        for (t = 0; t < 4; t++) {
            double acc = 0;
            for (j = -r; j <= r; j++) {
                int k = b + t * spb + j;
                if (k >= 0 && k < s->nb) acc += col[k];
            }
            pw[n][t] = (float)acc;
        }
    }
}

/* ------------------------------------------------------------------ */
/* soft metrics and decoding                                            */

typedef struct {
    float d[PI4_NSYM];       /* P(data=1 tone) - P(data=0 tone), local noise units */
    float p0[PI4_NSYM], p1[PI4_NSYM]; /* the two possible tones' powers       */
    double rho, sd;          /* signal excess, noise sd of d                 */
    int bessel;              /* 1: single-bin powers, exact noncoherent LLR   */
} soft;

/* ln I0(x), Abramowitz & Stegun 9.8.1 / 9.8.2 */
static double lnI0(double x)
{
    if (x < 0) x = -x;
    if (x <= 3.75) {
        double t = (x / 3.75) * (x / 3.75);
        return log(1 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813))))));
    } else {
        double t = 3.75 / x;
        return x - 0.5 * log(x) + log(0.39894228 + t * (0.01328592 + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 +
               t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377))))))));
    }
}

static void make_soft(float (*pw)[4], int bessel, soft *sf)
{
    int n; double rho = 0, v = 0, nl = 0;
    for (n = 0; n < PI4_NSYM; n++) nl += pw[n][1 - pi4_sync[n]] + pw[n][3 - pi4_sync[n]];
    nl /= 2 * PI4_NSYM;                      /* local noise level (impossible tones) */
    if (nl <= 0) nl = 1e-9;
    sf->bessel = bessel;
    for (n = 0; n < PI4_NSYM; n++) {
        int s = pi4_sync[n];
        float a = pw[n][s] / (float)nl, b = pw[n][s + 2] / (float)nl, c = pw[n][1 - s] / (float)nl, e = pw[n][3 - s] / (float)nl;
        sf->p0[n] = a; sf->p1[n] = b;
        sf->d[n] = b - a;
        rho += (a + b) - (c + e);
        v += (double)(c - e) * (c - e);
    }
    sf->rho = rho / PI4_NSYM;
    sf->sd = sqrt(v / PI4_NSYM);
    if (sf->sd <= 0) sf->sd = 1e-6;
    if (sf->rho < 0.05 * sf->sd) sf->rho = 0.05 * sf->sd;
}

#define MET_SCALE 32
#define MET_FLOOR (-20 * MET_SCALE)

static int16_t bit_metric(double l)
{
    /* log2(2/(1+e^-l)) - 1/2, l = LLR in favour of this bit value */
    double m;
    if (l < -30) m = l / 0.693147 - 0.5;
    else m = 1.0 - log(1 + exp(-l)) / 0.693147 - 0.5;
    m *= MET_SCALE;
    if (m < MET_FLOOR) m = MET_FLOOR;
    return (int16_t)floor(m + 0.5);
}

/* Fill Fano metrics from soft values with LLR scale g; returns nothing. */
/* Per-bit Fano metrics from one of three LLR models (powers in local noise
 * units, rho = estimated mean signal excess):
 *  MET_RAYLEIGH  L = (P1-P0) * rho/(1+rho): signal amplitude unknown per symbol
 *                (Rayleigh prior) - robust to fading, AGC, dropouts; the default
 *  MET_BESSEL    L = lnI0(2 sqrt(rho P1)) - lnI0(2 sqrt(rho P0)): exact for a
 *                constant-amplitude tone in AWGN, single-bin powers only
 *  MET_GAUSS     L = 2 rho (P1-P0) / var: Gaussian approximation, used for
 *                multi-bin (spread) windows */
enum { MET_RAYLEIGH, MET_BESSEL, MET_GAUSS };

static void build_metrics(const soft *sf, double scale, int kind, int16_t met[PI4_NSYM][2])
{
    int i;
    double rho = scale * sf->rho, g;
    if (kind == MET_BESSEL && !sf->bessel) kind = MET_GAUSS;
    g = kind == MET_RAYLEIGH ? rho / (1 + rho) * (sf->bessel ? 1.0 : 2.0 / (sf->sd * sf->sd))
                             : scale * 2 * sf->rho / (sf->sd * sf->sd);
    for (i = 0; i < PI4_NSYM; i++) {
        int n = pi4_interleave[i];
        double l = kind == MET_BESSEL ? lnI0(2 * sqrt(rho * sf->p1[n])) - lnI0(2 * sqrt(rho * sf->p0[n]))
                                      : g * sf->d[n];
        if (l > 12) l = 12;
        if (l < -12) l = -12;
        met[i][1] = bit_metric(l);
        met[i][0] = bit_metric(-l);
    }
}

typedef struct {
    uint64_t n; int ok; long metric; double fit, fit_z; int hard;
} attempt;

static void evaluate(const soft *sf, uint64_t n, attempt *a)
{
    uint8_t tones[PI4_NSYM]; int i; double num = 0, den = 0, ss = 0; int hard = 0;
    pi4_encode_n(n, tones);
    for (i = 0; i < PI4_NSYM; i++) {
        int bit = tones[i] >> 1;
        double sgn = bit ? 1 : -1;
        num += sgn * sf->d[i]; den += fabs(sf->d[i]); ss += (double)sf->d[i] * sf->d[i];
        if ((sf->d[i] > 0) != bit) hard++;
    }
    a->n = n; a->hard = hard;
    a->fit = den > 0 ? num / den : 0;
    a->fit_z = num / (sf->sd * sqrt((double)PI4_NSYM));
}

static int try_decode(const soft *sf, const pi4_params *p, double z, attempt *best)
{
    /* effort follows the sync statistic: a cell at z < 5.5 is almost never a
     * decodable signal (a -24 dB one scores 7-9), so it gets a quarter of the
     * cycle budget and the two unit-scale attempts only */
    long cycles_budget = z < 5.5 ? p->fano_cycles / 4 : p->fano_cycles;
    int nattempt = z < 5.5 ? 2 : 7;
    /* attempts: LLR scale x direction; the reversed code (bidirectional
     * sequential decoding) fails on different noise bursts than the forward */
    static const double scales[] = { 1.0, 1.0, 1.0, 0.5, 2.0, 0.5, 0.25 };
    static const int revs[]      = { 0,   1,   0,   0,   0,   1,   0 };
    static const int kinds[]     = { MET_RAYLEIGH, MET_RAYLEIGH, MET_BESSEL, MET_RAYLEIGH, MET_RAYLEIGH, MET_BESSEL, MET_RAYLEIGH };
    int k, found = 0;
    int16_t met[PI4_NSYM][2];
    best->ok = 0; best->fit_z = -1e9;
    for (k = 0; k < nattempt; k++) {
        uint8_t bits[PI4_NBITS];
        long metric = 0, cycles = 0;
        attempt a;
        int i;
        if (k == 0 || scales[k] != scales[k - 1] || kinds[k] != kinds[k - 1]) build_metrics(sf, scales[k], p->metric >= 0 ? p->metric : kinds[k], met);
        if (!pi4_fano(&met[0][0], 2 * MET_SCALE, cycles_budget, revs[k], bits, &metric, &cycles)) continue;
        a.n = 0;
        for (i = 0; i < PI4_NINFO; i++) a.n = (a.n << 1) | bits[i];
        evaluate(sf, a.n, &a);
        a.ok = 1; a.metric = metric;
        if (p->verbose)
            fprintf(stderr, "  fano %s scale %.2f%s: n=%llu fit=%.3f fit_z=%.1f hard=%d metric=%ld cycles=%ld\n",
                    kinds[k] == MET_RAYLEIGH ? "rayleigh" : kinds[k] == MET_BESSEL ? "bessel" : "gauss", scales[k], revs[k] ? " rev" : "", (unsigned long long)a.n, a.fit, a.fit_z, a.hard, metric, cycles);
        if (!found || a.fit_z > best->fit_z) { *best = a; found = 1; }
        if (a.fit_z > 12) break;      /* unambiguous */
    }
    return found;
}

/* ------------------------------------------------------------------ */

static double snr_of(float (*pw)[4], const uint8_t *tones, int navg, int wbins, double *energy)
{
    /* noise from the two tones each symbol's sync bit rules out (unbiased by
     * the signal's own leakage, unlike the whitening floor near a strong
     * tone); signal = decoded tone minus that noise */
    int n, t; double sig = 0, noise = 0;
    for (n = 0; n < PI4_NSYM; n++) {
        double other = 0;
        for (t = 0; t < 4; t++) if (t != tones[n]) other += pw[n][t];
        sig += pw[n][tones[n]] - other / 3;
        noise += pw[n][1 - pi4_sync[n]] + pw[n][3 - pi4_sync[n]];
    }
    *energy = sig;
    noise /= 2 * PI4_NSYM;
    sig /= PI4_NSYM;
    if (noise <= 0) noise = 1e-9;
    sig /= noise / wbins;                     /* over the per-bin noise density; both sides
                                                 scale with navg so the ratio is per minute */
    (void)navg;
    if (sig < 1e-3) sig = 1e-3;
    return 10 * log10(sig * 4.8 / PI4_SYMLEN);
}

static int overlaps(const pi4_result *a, const pi4_result *b)
{
    double sa = pi4_spacing_hz(a->variant), sb = pi4_spacing_hz(b->variant);
    double alo = a->tone0_hz - sa, ahi = a->tone0_hz + 4 * sa;
    double blo = b->tone0_hz - sb, bhi = b->tone0_hz + 4 * sb;
    if (fabs(a->t_s - b->t_s) > PI4_NSYM * PI4_SYMBOL_S) return 0;
    return alo < bhi && blo < ahi;
}

static int res_cmp(const void *a, const void *b)
{
    /* rank by explained signal energy (plausible texts first): an alias
     * hypothesis riding on one tone of a strong signal explains far less */
    const pi4_result *x = (const pi4_result *)a, *y = (const pi4_result *)b;
    double ex = x->energy / sqrt(x->spread_hz), ey = y->energy / sqrt(y->spread_hz);
    if (x->plausible != y->plausible) return y->plausible - x->plausible;
    return ex > ey ? -1 : ex < ey;
}

int pi4_decode_spec(const pi4_spec *s, const float *audio, size_t n, int fs,
                    const pi4_params *p, pi4_result *out, int max_out)
{
    cand cands[MAX_CAND]; int ncand = 0, i, v, w, nres = 0, tried = 0;
    float (*pw)[4] = (float (*)[4])malloc(sizeof(float) * 4 * PI4_NSYM);
    pi4_result *res = (pi4_result *)calloc((size_t)MAX_CAND, sizeof(pi4_result));
    float *x = NULL; size_t nx = 0;
    uint8_t (*ctones)[PI4_NSYM] = NULL;
    int nspread = p->max_spread + 1;
    clock_t t_start = clock();
    if (nspread > PI4_NSPREAD) nspread = PI4_NSPREAD;
    if (nspread < 1) nspread = 1;

    if (audio) {
        if (fs == PI4_FS) { x = (float *)audio; nx = n; }
        else x = pi4_resample12k(audio, n, fs, &nx);
    }
    if (p->ncalls > 0) {
        ctones = (uint8_t (*)[PI4_NSYM])malloc((size_t)p->ncalls * PI4_NSYM);
        for (i = 0; i < p->ncalls; i++)
            if (pi4_encode(p->calls[i], ctones[i])) memset(ctones[i], 0, PI4_NSYM);
    }

    /* 1. candidate search on three spread windows (1, 5, 17 bins); the
     *    exact window is chosen per candidate afterwards */
    for (w = 0; w < nspread; w += 2) {
        for (v = 0; v < 4; v++) {
            if (!(p->variants & (1u << v))) continue;
            if (pi4_spread_bins[w] >= 2 * pi4_variant_k((pi4_variant)v)) continue;
            search_cells(s, p, v, w, -1, NULL, p->min_sync_z, cands, &ncand);
            for (i = 0; i < p->ncalls; i++)
                search_cells(s, p, v, w, i, ctones[i], p->list_min_z, cands, &ncand);
        }
    }
    if (p->verbose) fprintf(stderr, "search: %d candidates, %.3f s\n", ncand, (double)(clock() - t_start) / CLOCKS_PER_SEC);
    qsort(cands, (size_t)ncand, sizeof(cand), cand_cmp);
    if (p->verbose)
        for (i = 0; i < ncand; i++)
            fprintf(stderr, "cand %2d: %s w=%d h=%d b=%d (%.1f Hz, t=%.2f s) z=%.1f call=%d\n", i,
                    pi4_variant_name((pi4_variant)cands[i].v), pi4_spread_bins[cands[i].w],
                    cands[i].h, cands[i].b, (s->b0 + cands[i].b + cands[i].db) * BIN_HZ,
                    (s->h0 + cands[i].h + cands[i].dh) * HOP / (double)PI4_FS, cands[i].z, cands[i].call);

    /* 2. refine + decode the best candidates */
    for (i = 0; i < ncand && tried < p->max_cand; i++) {
        cand *c = &cands[i];
        double sp = pi4_spacing_hz((pi4_variant)c->v);
        int spb = 2 * pi4_variant_k((pi4_variant)c->v);
        double f0 = (s->b0 + c->b + c->db) * BIN_HZ;
        double t0 = (s->h0 + c->h + c->dh) * HOP;      /* samples */
        soft sf; attempt a; pi4_result r; uint8_t tones[PI4_NSYM];
        int exact, wbest = c->w, k, d0 = 0, d1 = 0; double zbest;
        /* a cell inside the region of an already strong decode is that signal */
        for (k = 0; k < nres; k++) {
            double dsp = pi4_spacing_hz(res[k].variant);
            if (res[k].fit_z > 20 && fabs(res[k].t_s - t0 / PI4_FS) < PI4_NSYM * PI4_SYMBOL_S &&
                f0 - sp < res[k].tone0_hz + 4 * dsp && res[k].tone0_hz - dsp < f0 + 4 * sp) break;
        }
        if (k < nres) continue;
        tried++;
        zbest = -1e30;
        {
            /* At the edge of the search grid a wide window catches a signal just
             * outside the band with its skirt; left at the edge cell it would
             * report tone 0 up to half a window off, a false spread and a
             * noise-inflated S/N. There tone 0 is re-centred over every bin the
             * candidate's own window reached (blind and listed candidates alike:
             * sync_stat does not depend on the message). Inside the grid every
             * centre was already searched for that window: nothing moves. */
            int r0 = pi4_spread_bins[c->w] / 2, g0, g1;
            double fmin_, fmax_;
            band_for(p, &fmin_, &fmax_, &g0, &g1);
            g0 -= s->b0; g1 -= s->b0;                          /* the grid search_cells used */
            if (g0 < r0 + 1) g0 = r0 + 1;
            if (g1 > s->nb - 3 * spb - r0) g1 = s->nb - 3 * spb - r0;
            if (c->b - g0 < r0) d0 = -r0;
            if (g1 - 1 - c->b < r0) d1 = r0;
        }
        if (c->call < 0 || d0 != 0 || d1 != 0) {
            /* choose the spread window (largest self-normalised sync statistic),
             * at a grid edge jointly with tone 0 */
            int wi, db, bbest = c->b;
            for (db = d0; db <= d1; db++) {
                int bb = c->b + db;
                if (bb < 0 || bb + 3 * spb >= s->nb) continue;
                for (wi = 0; wi < nspread; wi++) {
                    double zz;
                    if (pi4_spread_bins[wi] >= spb) break;
                    spec_powers_w(s, c->h, bb, spb, pi4_spread_bins[wi], pw);
                    zz = sync_stat(pw);
                    if (zz > zbest) { zbest = zz; wbest = wi; bbest = bb; }
                }
            }
            /* the old cell's sub-bin offset means nothing at the new one */
            if (bbest != c->b) { c->b = bbest; c->db = 0; f0 = (s->b0 + c->b) * BIN_HZ; }
        }
        /* a narrow winner (1 or 3 bins) is refined on the time-domain audio at
         * the exact frequency, which beats any integer-bin window for an
         * unspread signal; the refined result is kept if it scores higher */
        exact = (x != NULL && wbest <= 1 && s->navg == 1);
        if (exact) {
            /* coordinate search: time (+-1/2 hop in 1/4 hops), then frequency (+-1/2 bin in 1/4 bins) */
            static const int tsteps[] = { -2, -1, 0, 1, 2 };
            double best = -1e30, bt = t0, bf = f0; int k;
            for (k = 0; k < 5; k++) {
                double tt = t0 + tsteps[k] * HOP / 4.0;
                double sc = exact_sync(s, x, nx, (long)floor(tt + 0.5), f0, sp, NULL);
                if (sc > best) { best = sc; bt = tt; }
            }
            for (k = 0; k < 5; k++) {
                double ff = f0 + tsteps[k] * BIN_HZ / 4.0;
                double sc = exact_sync(s, x, nx, (long)floor(bt + 0.5), ff, sp, NULL);
                if (sc > best) { best = sc; bf = ff; }
            }
            exact_sync(s, x, nx, (long)floor(bt + 0.5), bf, sp, pw);
            if (wbest == 0 || c->call >= 0 || sync_stat(pw) >= zbest) { t0 = bt; f0 = bf; wbest = 0; }
            else exact = 0;
        }
        if (!exact) spec_powers_w(s, c->h, c->b, spb, pi4_spread_bins[wbest], pw);
        c->w = wbest;
        make_soft(pw, wbest == 0 && s->navg == 1, &sf);
        if (p->verbose) fprintf(stderr, "  refined: t=%.3f s f0=%.2f Hz w=%d exact=%d rho=%.2f sd=%.2f\n", t0 / PI4_FS, f0, pi4_spread_bins[wbest], exact, sf.rho, sf.sd);
        if (p->verbose > 1) { uint8_t tt[PI4_NSYM]; int q; pi4_encode("OZ7IGY", tt);
            for (q = 0; q < PI4_NSYM; q++) fprintf(stderr, "    sym %3d tone %d: %8.1f %8.1f %8.1f %8.1f %s\n", q, tt[q], pw[q][0], pw[q][1], pw[q][2], pw[q][3],
                (pw[q][tt[q]] < pw[q][tt[q] ^ 2]) ? "ERR" : ""); }
        if (p->verbose) {
            int q; double mw = 0, vw = 0;
            for (q = 0; q < PI4_NSYM; q++) { double c1 = pw[q][1 - pi4_sync[q]], c2 = pw[q][3 - pi4_sync[q]]; mw += c1 + c2; vw += c1 * c1 + c2 * c2; }
            mw /= 2 * PI4_NSYM; vw = vw / (2 * PI4_NSYM) - mw * mw;
            fprintf(stderr, "  impossible tones: mean %.3f var %.3f (exact=%d w=%d)\n", mw, vw, exact, pi4_spread_bins[wbest]);
        }
        memset(&r, 0, sizeof r);
        if (c->call >= 0) {
            uint64_t nn; pi4_pack(p->calls[c->call], &nn);
            evaluate(&sf, nn, &a); a.ok = 1; a.metric = 0;
            r.listed = 1;
            if (a.fit_z < p->list_min_z * 0.7) continue;
        } else {
            if (!try_decode(&sf, p, c->z, &a)) continue;
        }
        pi4_unpack(a.n, r.text);
        r.plausible = pi4_text_plausible(r.text);
        r.variant = (pi4_variant)c->v;
        r.tone0_hz = f0; r.t_s = t0 / PI4_FS;
        r.edge = outside_band(p, r.variant, f0);
        r.sync_z = c->z; r.fit = a.fit; r.fit_z = a.fit_z; r.hard_errors = a.hard;
        r.spread_hz = pi4_spread_bins[c->w] * BIN_HZ;
        r.fano_metric = a.metric;
        pi4_encode_n(a.n, tones);
        r.snr_db = snr_of(pw, tones, s->navg, pi4_spread_bins[wbest], &r.energy);
        if (a.n == 0) continue;                   /* all-zero codeword: never a real message */
        if (p->verbose)
            fprintf(stderr, "  -> %-8s fit_z=%.1f fit=%.3f hard=%d snr=%.1f plausible=%d rho=%.2f sd=%.2f\n",
                    r.text, r.fit_z, r.fit, r.hard_errors, r.snr_db, r.plausible, sf.rho, sf.sd);
        if (!r.listed) {
            double gate = p->min_fit_z;
            if (!r.plausible) gate += 3.0;
            if (p->relax) gate -= 1.0;
            if (r.fit_z < gate) continue;
            if (r.hard_errors > 60) continue;
        }
        res[nres++] = r;
    }

    if (p->verbose) fprintf(stderr, "decode: %d tried, %d accepted, %.3f s\n", tried, nres, (double)(clock() - t_start) / CLOCKS_PER_SEC);
    /* 3. de-duplicate: one decode per overlapping time/frequency region */
    qsort(res, (size_t)nres, sizeof(pi4_result), res_cmp);
    {
        int k, m = 0;
        for (i = 0; i < nres && m < max_out; i++) {
            int dup = 0;
            for (k = 0; k < m; k++) if (overlaps(&out[k], &res[i])) { dup = 1; break; }
            if (!dup) out[m++] = res[i];
        }
        nres = m;
    }
    if (x && x != audio) free(x);
    free(pw); free(res); free(ctones);
    return nres;
}

int pi4_decode(const float *audio, size_t n, int fs, const pi4_params *p,
               pi4_result *out, int max_out)
{
    pi4_spec *s;
    int r;
    float *x = NULL; size_t nx = n;
    if (fs != PI4_FS) { x = pi4_resample12k(audio, n, fs, &nx); if (!x) return 0; audio = x; fs = PI4_FS; }
    s = pi4_spec_build(audio, nx, fs, p);
    if (!s) { free(x); return 0; }
    r = pi4_decode_spec(s, audio, nx, fs, p, out, max_out);
    pi4_spec_free(s);
    free(x);
    return r;
}
