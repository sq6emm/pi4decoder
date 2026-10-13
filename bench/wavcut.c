/* wavcut in.wav out.wav skip_s [len_s]: drop the first skip_s seconds (and
 * optionally truncate), for reference decoders that want the minute at the
 * start of the file. SPDX-License-Identifier: GPL-3.0-or-later */
#include "../src/wav.h"
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
    int fs; size_t n, a, b;
    float *x;
    if (argc < 4) { fprintf(stderr, "usage: wavcut in.wav out.wav skip_s [len_s]\n"); return 2; }
    x = wav_read(argv[1], &fs, &n);
    if (!x) { fprintf(stderr, "cannot read %s\n", argv[1]); return 1; }
    a = (size_t)(atof(argv[3]) * fs); if (a > n) a = n;
    b = argc > 4 ? a + (size_t)(atof(argv[4]) * fs) : n; if (b > n) b = n;
    return wav_write(argv[2], x + a, b - a, fs) ? 1 : 0;
}
