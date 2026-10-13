/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef WAV_H
#define WAV_H
#include <stddef.h>
#include <stdint.h>
float *wav_read(const char *path, int *fs, size_t *n);   /* malloc'd mono float */
int    wav_write(const char *path, const float *x, size_t n, int fs); /* 16-bit PCM */
#endif
