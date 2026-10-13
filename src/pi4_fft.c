/*
 * Minimal radix-2 FFT (complex, in place) and a real-input wrapper that
 * computes an N-point real FFT with an N/2-point complex FFT. Enough for a
 * 4096-point spectrogram column; no external library.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "pi4_fft.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int pi4_fft_init(pi4_fft *f, int n)
{
    int i, j, k;
    f->n = n;
    f->tw = (float *)malloc(sizeof(float) * (size_t)n);      /* cos,sin pairs for n/2 */
    f->tw2 = (float *)malloc(sizeof(float) * (size_t)n * 2);  /* for the real split */
    f->rev = (uint16_t *)malloc(sizeof(uint16_t) * (size_t)n);
    if (!f->tw || !f->tw2 || !f->rev) return -1;
    for (i = 0; i < n / 2; i++) {
        f->tw[2 * i] = (float)cos(2 * M_PI * i / n);
        f->tw[2 * i + 1] = (float)-sin(2 * M_PI * i / n);
    }
    for (i = 0; i < n; i++) {
        f->tw2[2 * i] = (float)cos(M_PI * i / n);
        f->tw2[2 * i + 1] = (float)-sin(M_PI * i / n);
    }
    for (i = 0; i < n; i++) {
        for (j = 0, k = i; (1 << ++j) < n; ) ;
        k = 0;
        for (j = 0; (1 << j) < n; j++) k = (k << 1) | ((i >> j) & 1);
        f->rev[i] = (uint16_t)k;
    }
    return 0;
}

void pi4_fft_free(pi4_fft *f)
{
    free(f->tw); free(f->tw2); free(f->rev);
    f->tw = f->tw2 = NULL; f->rev = NULL;
}

/* In-place complex FFT of x[2*n] (re,im interleaved). */
void pi4_fft_complex(const pi4_fft *f, float *x)
{
    int n = f->n, i, len;
    for (i = 0; i < n; i++) {
        int j = f->rev[i];
        if (j > i) {
            float tr = x[2 * i], ti = x[2 * i + 1];
            x[2 * i] = x[2 * j]; x[2 * i + 1] = x[2 * j + 1];
            x[2 * j] = tr; x[2 * j + 1] = ti;
        }
    }
    for (len = 2; len <= n; len <<= 1) {
        int half = len >> 1, step = n / len, j;
        for (i = 0; i < n; i += len) {
            for (j = 0; j < half; j++) {
                const float wr = f->tw[2 * j * step], wi = f->tw[2 * j * step + 1];
                float *a = x + 2 * (i + j), *b = x + 2 * (i + j + half);
                float br = b[0] * wr - b[1] * wi, bi = b[0] * wi + b[1] * wr;
                b[0] = a[0] - br; b[1] = a[1] - bi;
                a[0] += br; a[1] += bi;
            }
        }
    }
}

/* Real FFT of 2n real samples in[0..2n-1] (in is packed as n complex
 * values, modified) -> out[k] = X[k] for k = 0..n (re,im pairs, n+1 bins). */
void pi4_fft_real(const pi4_fft *f, float *in, float *out)
{
    int n = f->n, k;
    pi4_fft_complex(f, in);
    for (k = 0; k <= n; k++) {
        int k1 = k & (n - 1), k2 = (n - k) & (n - 1);
        float ar = in[2 * k1], ai = in[2 * k1 + 1];
        float br = in[2 * k2], bi = -in[2 * k2 + 1];     /* conj(Z[n-k]) */
        float er = 0.5f * (ar + br), ei = 0.5f * (ai + bi);     /* even part */
        float or_ = 0.5f * (ai - bi), oi = -0.5f * (ar - br);   /* odd/(j)   */
        float wr = f->tw2[2 * k], wi = f->tw2[2 * k + 1];
        out[2 * k] = er + or_ * wr - oi * wi;
        out[2 * k + 1] = ei + or_ * wi + oi * wr;
    }
}
