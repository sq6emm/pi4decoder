/* False-decode measurement on pure Gaussian noise minutes.
 *
 *   test_noise [minutes] [seed] [gate] [extra: --calls LIST | --all-variants]
 *
 * Runs the blind decoder with the given fit_z gate (default: library default)
 * on `minutes` independent noise minutes (all four variants, +-150 Hz, the
 * default monitor configuration) and prints every decode that passes the
 * gate plus the maximum fit_z / sync_z seen with the gate disabled, so the
 * gate margin can be read off. Exit status 1 if anything passed the gate.
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "pi4.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static uint64_t rs;
static double uniform(void)
{
    rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27;
    return (((rs * 0x2545F4914F6CDD1Dull) >> 32) + 0.5) / 4294967296.0;
}
static double gauss(void) { return sqrt(-2 * log(uniform())) * cos(2 * M_PI * uniform()); }

int main(int argc, char **argv)
{
    int minutes = argc > 1 ? atoi(argv[1]) : 100;
    unsigned long seed = argc > 2 ? strtoul(argv[2], NULL, 10) : 1;
    double gate = argc > 3 ? atof(argv[3]) : -1;
    const char *calls[32]; int ncalls = 0; char *callbuf = NULL;
    int i, passed = 0, converged = 0;
    double max_fit_z = 0, max_sync_z = 0, max_fit = 0;
    size_t n = 12000 * 60, k;
    float *x = (float *)malloc(sizeof(float) * n);
    pi4_params p, p0;
    clock_t t0 = clock();

    for (i = 4; i < argc; i++)
        if (!strcmp(argv[i], "--calls") && i + 1 < argc) {
            char *t; callbuf = argv[++i];
            for (t = strtok(callbuf, ","); t && ncalls < 32; t = strtok(NULL, ",")) calls[ncalls++] = t;
        }
    pi4_params_default(&p);
    if (gate >= 0) p.min_fit_z = gate;
    p.t_lo_s = 0.5; p.t_hi_s = 5.5;
    p.calls = calls; p.ncalls = ncalls;
    p0 = p; p0.min_fit_z = 0; p0.relax = 1;      /* gate off: see everything that converges */
    for (i = 0; i < minutes; i++) {
        pi4_result r[8]; int nr, j;
        rs = 0x9E3779B97F4A7C15ull ^ ((uint64_t)seed * 1000003u + (uint64_t)i * 7919u);
        for (k = 0; k < 100; k++) uniform();
        for (k = 0; k < n; k++) x[k] = (float)gauss();
        nr = pi4_decode(x, n, 12000, &p0, r, 8);
        for (j = 0; j < nr; j++) {
            converged++;
            if (r[j].fit_z > max_fit_z) max_fit_z = r[j].fit_z;
            if (r[j].sync_z > max_sync_z) max_sync_z = r[j].sync_z;
            if (r[j].fit > max_fit) max_fit = r[j].fit;
            if (r[j].listed || r[j].fit_z >= p.min_fit_z + (r[j].plausible ? 0 : 3.0)) {
                passed++;
                printf("FALSE DECODE minute %d: '%s' %s tone0=%.1f fit_z=%.1f fit=%.3f sync_z=%.1f hard=%d plausible=%d listed=%d\n",
                       i, r[j].text, pi4_variant_name(r[j].variant), r[j].tone0_hz, r[j].fit_z, r[j].fit,
                       r[j].sync_z, r[j].hard_errors, r[j].plausible, r[j].listed);
            } else {
                printf("  converged (rejected) minute %d: '%s' fit_z=%.1f fit=%.3f sync_z=%.1f hard=%d plausible=%d\n",
                       i, r[j].text, r[j].fit_z, r[j].fit, r[j].sync_z, r[j].hard_errors, r[j].plausible);
            }
        }
    }
    printf("test_noise: %d minutes, seed %lu, gate %.1f: %d converged, %d passed; max fit_z %.2f, max fit %.3f, max sync_z %.2f; %.2f s/minute\n",
           minutes, seed, p.min_fit_z, converged, passed, max_fit_z, max_fit, max_sync_z,
           (double)(clock() - t0) / CLOCKS_PER_SEC / minutes);
    free(x);
    return passed != 0;
}
