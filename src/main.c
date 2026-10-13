/*
 * pi4decoder command line.
 *
 *   pi4decoder [options] file.wav [boundary_s]
 *
 * Prints one JSON object per decode. See README.md and --help.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "pi4.h"
#include "wav.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(void)
{
    fprintf(stderr,
        "usage: pi4decoder [options] file.wav [boundary_s] [file2.wav boundary2_s ...]\n"
        "  boundary_s      where the even UTC minute (PI4 start) falls in the file (default 0)\n"
        "  --radius HZ     tone-0 search +-HZ around the convention (default 150)\n"
        "  --center HZ     CW carrier convention: tone0 = center - spacing/2 (default 800)\n"
        "  --band LO HI    search tone 0 in LO..HI Hz instead of center/radius\n"
        "  --variants L    comma list of PI4,PI4-80,PI4-96,PI4-120 (default all)\n"
        "  --dt S          message start searched +-S s around the boundary (default 2.5)\n"
        "  --scan          search the whole file for any message start (boundary ignored;\n"
        "                  default band 300..2700 Hz unless --band/--radius given)\n"
        "  --spread N      widest tone energy window: 0=2.9 Hz .. 5=97 Hz (default 5)\n"
        "  --calls LIST    also match these known messages (comma list), reported with listed=1\n"
        "  --avg           average all given files (same alignment) before decoding\n"
        "  --gate Z        blind-decode acceptance gate on fit_z (default 6.0)\n"
        "  --relax         accept marginal decodes one unit below the gate\n"
        "  --cycles N      Fano cycle budget per attempt (default 200000)\n"
        "  --max-cand N    candidate cells tried (default 12)\n"
        "  --min-sync Z    candidate gate on the sync z-score (default 4.0)\n"
        "  --fs HZ         override the WAV sample rate (clock-error experiments)\n"
        "  --metric K      force an LLR model: 0 rayleigh, 1 bessel, 2 gauss (testing)\n"
        "  --verbose       diagnostics on stderr (twice: per-symbol dump)\n");
}

static void json_str(const char *s)
{
    putchar('"');
    for (; *s; s++) {
        if (*s == '"' || *s == '\\') putchar('\\');
        putchar(*s);
    }
    putchar('"');
}

static void print_result(const pi4_result *r, double boundary_s)
{
    printf("{\"text\": "); json_str(r->text);
    printf(", \"variant\": \"%s\", \"tone0_hz\": %.2f, \"dt_s\": %.3f, \"t_s\": %.3f, \"snr_db\": %.1f, "
           "\"fit\": %.3f, \"fit_z\": %.1f, \"sync_z\": %.1f, \"hard_errors\": %d, \"spread_hz\": %.1f, "
           "\"plausible\": %s, \"listed\": %s}\n",
           pi4_variant_name(r->variant), r->tone0_hz, r->t_s - boundary_s, r->t_s, r->snr_db,
           r->fit, r->fit_z, r->sync_z, r->hard_errors, r->spread_hz,
           r->plausible ? "true" : "false", r->listed ? "true" : "false");
}

int main(int argc, char **argv)
{
    pi4_params p;
    const char *files[64]; double bounds[64]; int nfiles = 0;
    const char *calls[64]; int ncalls = 0; char *callbuf = NULL;
    int scan = 0, avg = 0, i, band_set = 0, radius_set = 0, fs_override = 0, maxcand_set = 0;
    double dt = 2.5;
    pi4_result res[16];

    pi4_params_default(&p);
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--radius") && i + 1 < argc) { p.radius_hz = atof(argv[++i]); radius_set = 1; }
        else if (!strcmp(a, "--center") && i + 1 < argc) p.center_hz = atof(argv[++i]);
        else if (!strcmp(a, "--band") && i + 2 < argc) { p.fmin_hz = atof(argv[++i]); p.fmax_hz = atof(argv[++i]); band_set = 1; }
        else if (!strcmp(a, "--variants") && i + 1 < argc) {
            char *s = argv[++i], *t;
            p.variants = 0;
            for (t = strtok(s, ","); t; t = strtok(NULL, ",")) {
                int v = pi4_variant_parse(t);
                if (v < 0) { fprintf(stderr, "unknown variant %s\n", t); return 2; }
                p.variants |= 1u << v;
            }
        }
        else if (!strcmp(a, "--dt") && i + 1 < argc) dt = atof(argv[++i]);
        else if (!strcmp(a, "--scan")) scan = 1;
        else if (!strcmp(a, "--avg")) avg = 1;
        else if (!strcmp(a, "--spread") && i + 1 < argc) p.max_spread = atoi(argv[++i]);
        else if (!strcmp(a, "--gate") && i + 1 < argc) p.min_fit_z = atof(argv[++i]);
        else if (!strcmp(a, "--relax")) p.relax = 1;
        else if (!strcmp(a, "--cycles") && i + 1 < argc) p.fano_cycles = atol(argv[++i]);
        else if (!strcmp(a, "--max-cand") && i + 1 < argc) { p.max_cand = atoi(argv[++i]); maxcand_set = 1; }
        else if (!strcmp(a, "--min-sync") && i + 1 < argc) p.min_sync_z = atof(argv[++i]);
        else if (!strcmp(a, "--list-gate") && i + 1 < argc) p.list_min_z = atof(argv[++i]);
        else if (!strcmp(a, "--calls") && i + 1 < argc) {
            char *t;
            callbuf = argv[++i];
            for (t = strtok(callbuf, ","); t && ncalls < 64; t = strtok(NULL, ",")) calls[ncalls++] = t;
        }
        else if (!strcmp(a, "--fs") && i + 1 < argc) fs_override = atoi(argv[++i]);
        else if (!strcmp(a, "--metric") && i + 1 < argc) p.metric = atoi(argv[++i]);
        else if (!strcmp(a, "--verbose") || !strcmp(a, "-v")) p.verbose++;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); return 0; }
        else if (a[0] == '-' && a[1] == '-') { fprintf(stderr, "unknown option %s\n", a); usage(); return 2; }
        else if (nfiles < 64) {
            char *end;
            double b = strtod(a, &end);
            if (*end == 0 && nfiles > 0 && bounds[nfiles - 1] == -1) bounds[nfiles - 1] = b;
            else { files[nfiles] = a; bounds[nfiles] = -1; nfiles++; }
        }
    }
    if (nfiles == 0) { usage(); return 2; }
    for (i = 0; i < nfiles; i++) if (bounds[i] < 0) bounds[i] = 0;
    p.calls = calls; p.ncalls = ncalls;
    if (scan && !band_set && !radius_set) { p.fmin_hz = 300; p.fmax_hz = 2700; }
    if (scan && !maxcand_set) p.max_cand = 24;

    if (avg) {
        pi4_spec *acc = NULL;
        for (i = 0; i < nfiles; i++) {
            int fs; size_t n; float *x = wav_read(files[i], &fs, &n); pi4_spec *s;
            if (!x) { fprintf(stderr, "cannot read %s\n", files[i]); return 1; }
            p.t_lo_s = bounds[i] - dt; p.t_hi_s = bounds[i] + dt;
            s = pi4_spec_build(x, n, fs, &p);
            free(x);
            if (!acc) acc = s;
            else { if (pi4_spec_add(acc, s)) { fprintf(stderr, "%s: alignment differs, skipped\n", files[i]); } pi4_spec_free(s); }
        }
        p.t_lo_s = bounds[0] - dt; p.t_hi_s = bounds[0] + dt;
        {
            int n = pi4_decode_spec(acc, NULL, 0, 0, &p, res, 16), k;
            for (k = 0; k < n; k++) print_result(&res[k], bounds[0]);
        }
        pi4_spec_free(acc);
        return 0;
    }
    for (i = 0; i < nfiles; i++) {
        int fs; size_t n, k; float *x = wav_read(files[i], &fs, &n); int nr;
        if (!x) { fprintf(stderr, "cannot read %s\n", files[i]); return 1; }
        if (fs_override) fs = fs_override;
        if (scan) { p.t_lo_s = 0; p.t_hi_s = (double)n / fs - 20.0; if (p.t_hi_s < 0) p.t_hi_s = 0; }
        else { p.t_lo_s = bounds[i] - dt; p.t_hi_s = bounds[i] + dt; }
        nr = pi4_decode(x, n, fs, &p, res, 16);
        for (k = 0; k < (size_t)nr; k++) print_result(&res[k], scan ? 0 : bounds[i]);
        free(x);
    }
    return 0;
}
