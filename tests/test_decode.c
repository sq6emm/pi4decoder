/* Decoder tests: FFT correctness, resampler, synthetic end-to-end decodes
 * at several sample rates, noise produces nothing.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "pi4.h"
#include "pi4_fft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int fails = 0;
#define CHECK(c, ...) do { if (!(c)) { fails++; printf("FAIL %s:%d ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static uint64_t rs = 0x1234567887654321ull;
static double uniform(void)
{
    rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27;
    return (((rs * 0x2545F4914F6CDD1Dull) >> 32) + 0.5) / 4294967296.0;
}
static double gauss(void) { return sqrt(-2 * log(uniform())) * cos(2 * M_PI * uniform()); }

/* synthetic PI4 at fs, message start at lead s, S/N (2500 Hz) snr dB */
static float *synth(const char *text, int fs, double lead, double f0, double snr, double dur, size_t *n)
{
    uint8_t tones[PI4_NSYM];
    double sp = pi4_spacing_hz(PI4_V_PI4), ph = 0, symlen = fs / 6.0;
    double sd = sqrt(0.5 / pow(10, snr / 10) / (2500.0 / (fs / 2.0)));
    size_t k, N = (size_t)(dur * fs);
    float *x = (float *)calloc(N, sizeof(float));
    pi4_encode(text, tones);
    for (k = 0; k < N; k++) {
        double t = (double)k - lead * fs, s = 0;
        if (t >= 0 && t < PI4_NSYM * symlen) {
            double f = f0 + tones[(int)(t / symlen)] * sp;
            ph += 2 * M_PI * f / fs;
            s = cos(ph);
        }
        x[k] = (float)(s + sd * gauss());
    }
    *n = N;
    return x;
}

int main(void)
{
    /* 1. FFT against a naive DFT */
    {
        pi4_fft f; int N = 2048, k, n;
        float *in = (float *)malloc(sizeof(float) * (size_t)N), *out = (float *)malloc(sizeof(float) * (size_t)(N + 2));
        double maxerr = 0, maxmag = 0;
        pi4_fft_init(&f, N / 2);
        for (n = 0; n < N; n++) in[n] = (float)(gauss() + (n < 1000 ? cos(2 * M_PI * 100.3 * n / N) : 0));
        {
            float *copy = (float *)malloc(sizeof(float) * (size_t)N);
            memcpy(copy, in, sizeof(float) * (size_t)N);
            pi4_fft_real(&f, copy, out);
            for (k = 0; k <= N / 2; k += 7) {
                double re = 0, im = 0;
                for (n = 0; n < N; n++) { re += in[n] * cos(2 * M_PI * k * n / N); im -= in[n] * sin(2 * M_PI * k * n / N); }
                if (fabs(re - out[2 * k]) > maxerr) maxerr = fabs(re - out[2 * k]);
                if (fabs(im - out[2 * k + 1]) > maxerr) maxerr = fabs(im - out[2 * k + 1]);
                if (fabs(re) > maxmag) maxmag = fabs(re);
            }
            free(copy);
        }
        CHECK(maxerr < 1e-2 * maxmag, "real FFT vs DFT: max err %g (max %g)", maxerr, maxmag);
        pi4_fft_free(&f); free(in); free(out);
    }

    /* 2. resampler keeps a tone's frequency and amplitude */
    {
        int fs_in[] = { 6000, 8000, 24000, 48000 }, i;
        for (i = 0; i < 4; i++) {
            size_t n = (size_t)fs_in[i] * 2, m, k;
            float *x = (float *)malloc(sizeof(float) * n), *y;
            double re = 0, im = 0, amp;
            for (k = 0; k < n; k++) x[k] = (float)cos(2 * M_PI * 1000.0 * k / fs_in[i]);
            y = pi4_resample12k(x, n, fs_in[i], &m);
            CHECK(m >= 23900 && m <= 24000, "resampled length %zu from %d", m, fs_in[i]);
            for (k = 2000; k < 22000; k++) { re += y[k] * cos(2 * M_PI * 1000.0 * k / 12000); im -= y[k] * sin(2 * M_PI * 1000.0 * k / 12000); }
            amp = 2 * sqrt(re * re + im * im) / 20000;
            CHECK(fabs(amp - 1.0) < 0.05, "resample %d: tone amplitude %.3f", fs_in[i], amp);
            free(x); free(y);
        }
    }

    /* 3. end-to-end decodes at several rates and offsets */
    {
        struct { int fs; double lead, f0, snr; const char *text; } cases[] = {
            { 12000, 3.0, 682.8125, -10, "OZ7IGY" },
            { 12000, 3.7, 690.0, -18, "SR6LEG" },
            { 12000, 1.4, 650.0, -14, "ED5YAE/B" },
            { 6000, 3.0, 700.0, -12, "OZ7IGY" },
            { 8000, 2.2, 682.8125, -12, "IW9GDC/B" },
            { 24000, 3.0, 682.8125, -16, "G3ZJO" },
            { 48000, 3.0, 682.8125, -10, "PE1ITR/B" },
        };
        unsigned i;
        for (i = 0; i < sizeof cases / sizeof cases[0]; i++) {
            size_t n; float *x = synth(cases[i].text, cases[i].fs, cases[i].lead, cases[i].f0, cases[i].snr, 40, &n);
            pi4_params p; pi4_result r[4]; int nr;
            pi4_params_default(&p);
            p.t_lo_s = 3.0 - 2.5; p.t_hi_s = 3.0 + 2.5;
            nr = pi4_decode(x, n, cases[i].fs, &p, r, 4);
            CHECK(nr == 1 && !strcmp(r[0].text, cases[i].text), "case %u (%d Hz, %s, %g dB): %d decodes%s%s", i,
                  cases[i].fs, cases[i].text, cases[i].snr, nr, nr ? ": " : "", nr ? r[0].text : "");
            if (nr == 1) {
                CHECK(fabs(r[0].t_s - cases[i].lead) < 0.03, "case %u: t_s %.3f", i, r[0].t_s);
                CHECK(fabs(r[0].tone0_hz - cases[i].f0) < 1.5, "case %u: tone0 %.2f", i, r[0].tone0_hz);
                CHECK(fabs(r[0].snr_db - cases[i].snr) < 2.0, "case %u: snr %.1f vs %g", i, r[0].snr_db, cases[i].snr);
                CHECK(r[0].variant == PI4_V_PI4, "case %u: variant", i);
            }
            free(x);
        }
    }

    /* 4. pure noise minutes decode to nothing (all variants, full radius) */
    {
        int i, bad = 0;
        for (i = 0; i < 6; i++) {
            size_t n = 12000 * 60, k; float *x = (float *)malloc(sizeof(float) * n);
            pi4_params p; pi4_result r[4]; int nr;
            for (k = 0; k < n; k++) x[k] = (float)gauss();
            pi4_params_default(&p);
            p.t_lo_s = 0.5; p.t_hi_s = 5.5;
            nr = pi4_decode(x, n, 12000, &p, r, 4);
            if (nr) { bad++; printf("  noise %d decoded as '%s' fit_z=%.1f\n", i, r[0].text, r[0].fit_z); }
            free(x);
        }
        CHECK(bad == 0, "%d of 6 noise minutes produced a decode", bad);
    }

    printf(fails ? "test_decode: %d FAILURES\n" : "test_decode: all passed\n", fails);
    return fails != 0;
}
