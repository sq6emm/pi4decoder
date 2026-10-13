/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef PI4_FFT_H
#define PI4_FFT_H
#include <stdint.h>

typedef struct {
    int n;          /* complex size (real FFT length is 2n) */
    float *tw;      /* n/2 twiddles (re,im) */
    float *tw2;     /* n twiddles for the real split */
    uint16_t *rev;  /* bit reversal table */
} pi4_fft;

int  pi4_fft_init(pi4_fft *f, int n);
void pi4_fft_free(pi4_fft *f);
void pi4_fft_complex(const pi4_fft *f, float *x);
void pi4_fft_real(const pi4_fft *f, float *in, float *out);

#endif
