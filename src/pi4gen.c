/*
 * pi4gen: synthetic PI4 test signals (and pure noise) as 16-bit WAV.
 *
 *   pi4gen [options] out.wav
 *     --text MSG       message (default OZ7IGY)
 *     --snr DB         S/N in 2500 Hz (default 0 = clean); use --noise for noise only
 *     --fs HZ          sample rate (default 12000)
 *     --f0 HZ          tone-0 frequency (default 682.8125 = 800 - 117.1875)
 *     --variant V      PI4, PI4-80, PI4-96, PI4-120
 *     --lead S         message start time in the file (default 3.0)
 *     --dur S          file length (default 63)
 *     --spread HZ      Doppler spread (Rayleigh fading, Gaussian spectrum of this width)
 *     --drift HZ       linear frequency drift over the message
 *     --spacing F      tone-spacing error factor (1.017 for the mp4 case)
 *     --carrier        add the CW-ID/carrier at f0+spacing/2 after the message
 *     --noise          pure noise (no signal)
 *     --seed N
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "pi4.h"
#include "wav.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint64_t rng_state = 0x9E3779B97F4A7C15ull;
static uint32_t rng32(void)
{
    rng_state ^= rng_state >> 12; rng_state ^= rng_state << 25; rng_state ^= rng_state >> 27;
    return (uint32_t)((rng_state * 0x2545F4914F6CDD1Dull) >> 32);
}
static double uniform(void) { return (rng32() + 0.5) / 4294967296.0; }
static double gauss(void)
{
    double u1 = uniform(), u2 = uniform();
    return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}

int main(int argc, char **argv)
{
    const char *text = "OZ7IGY", *out = NULL;
    double snr = 0, f0 = 800 - 234.375 / 2, lead = 3.0, dur = 63, spread = 0, drift = 0, spfac = 1.0;
    int fs = 12000, variant = 0, carrier = 0, noise_only = 0, i, f0_set = 0;
    uint8_t tones[PI4_NSYM];
    size_t n, k;
    float *x;
    double sp, ph = 0;

    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--text") && i + 1 < argc) text = argv[++i];
        else if (!strcmp(a, "--snr") && i + 1 < argc) snr = atof(argv[++i]);
        else if (!strcmp(a, "--fs") && i + 1 < argc) fs = atoi(argv[++i]);
        else if (!strcmp(a, "--f0") && i + 1 < argc) { f0 = atof(argv[++i]); f0_set = 1; }
        else if (!strcmp(a, "--variant") && i + 1 < argc) { variant = pi4_variant_parse(argv[++i]); if (variant < 0) return 2; }
        else if (!strcmp(a, "--lead") && i + 1 < argc) lead = atof(argv[++i]);
        else if (!strcmp(a, "--dur") && i + 1 < argc) dur = atof(argv[++i]);
        else if (!strcmp(a, "--spread") && i + 1 < argc) spread = atof(argv[++i]);
        else if (!strcmp(a, "--drift") && i + 1 < argc) drift = atof(argv[++i]);
        else if (!strcmp(a, "--spacing") && i + 1 < argc) spfac = atof(argv[++i]);
        else if (!strcmp(a, "--carrier")) carrier = 1;
        else if (!strcmp(a, "--noise")) noise_only = 1;
        else if (!strcmp(a, "--seed") && i + 1 < argc) rng_state ^= (uint64_t)strtoull(argv[++i], NULL, 10) * 0x9E3779B97F4A7C15ull + 1;
        else if (a[0] != '-') out = a;
        else { fprintf(stderr, "unknown option %s\n", a); return 2; }
    }
    if (!out) { fprintf(stderr, "usage: pi4gen [options] out.wav\n"); return 2; }
    sp = pi4_spacing_hz((pi4_variant)variant) * spfac;
    if (!f0_set) f0 = 800 - pi4_spacing_hz((pi4_variant)variant) / 2;
    if (pi4_encode(text, tones)) { fprintf(stderr, "bad message '%s'\n", text); return 2; }
    for (i = 0; i < 1000; i++) rng32();

    n = (size_t)(dur * fs);
    x = (float *)calloc(n, sizeof(float));
    if (!noise_only) {
        /* fading process: sum of 24 complex sinusoids with Gaussian-distributed
         * frequencies of sd = spread/2 (so the spectrum's -3 dB width is ~spread) */
        const int NF = 24;
        double ff[24], fp[24];
        double symlen = fs / 6.0, t0 = lead * fs, msglen = PI4_NSYM * symlen;
        for (i = 0; i < NF; i++) { ff[i] = gauss() * spread / 2.355; fp[i] = 2 * M_PI * uniform(); }
        for (k = 0; k < n; k++) {
            double t = (double)k - t0;
            double f, gi, gq, s;
            int sym;
            if (t < 0) continue;
            if (t < msglen) {
                sym = (int)(t / symlen);
                f = f0 + tones[sym] * sp + drift * (t / msglen);
            } else if (carrier && t < 60.0 * fs) {
                f = f0 + sp / 2 + drift;
            } else break;
            ph += 2 * M_PI * f / fs;
            if (ph > 2 * M_PI) ph -= 2 * M_PI;
            if (spread > 0) {
                gi = 0; gq = 0;
                for (i = 0; i < NF; i++) {
                    double a = 2 * M_PI * ff[i] * k / fs + fp[i];
                    gi += cos(a); gq += sin(a);
                }
                gi /= sqrt((double)NF); gq /= sqrt((double)NF);
                s = gi * cos(ph) - gq * sin(ph);
            } else s = cos(ph);
            x[k] = (float)s;
        }
    }
    if (snr != 0 || noise_only) {
        double sig2 = 0.5 / pow(10, snr / 10) / (2500.0 / (fs / 2.0));
        double sd = sqrt(sig2);
        if (noise_only) sd = 1.0;
        for (k = 0; k < n; k++) x[k] += (float)(sd * gauss());
    }
    {   float peak = 0;
        for (k = 0; k < n; k++) if (fabsf(x[k]) > peak) peak = fabsf(x[k]);
        if (peak > 0) for (k = 0; k < n; k++) x[k] = x[k] / peak * 0.5f;
    }
    if (wav_write(out, x, n, fs)) { fprintf(stderr, "cannot write %s\n", out); return 1; }
    free(x);
    return 0;
}
