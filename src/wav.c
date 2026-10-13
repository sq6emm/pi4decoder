/* Minimal RIFF/WAVE reader and writer (16-bit PCM, mono or first channel).
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "wav.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t rd32(const unsigned char *p) { return p[0] | (p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); }
static uint16_t rd16(const unsigned char *p) { return (uint16_t)(p[0] | (p[1] << 8)); }

float *wav_read(const char *path, int *fs, size_t *n)
{
    FILE *fp = fopen(path, "rb");
    unsigned char hdr[12], ch[8];
    int channels = 1, bits = 16, fmt = 1;
    float *out = NULL;
    if (!fp) return NULL;
    if (fread(hdr, 1, 12, fp) != 12 || memcmp(hdr, "RIFF", 4) || memcmp(hdr + 8, "WAVE", 4)) { fclose(fp); return NULL; }
    *fs = 0;
    while (fread(ch, 1, 8, fp) == 8) {
        uint32_t len = rd32(ch + 4);
        if (!memcmp(ch, "fmt ", 4)) {
            unsigned char f[16];
            if (len < 16 || fread(f, 1, 16, fp) != 16) break;
            fmt = rd16(f); channels = rd16(f + 2); *fs = (int)rd32(f + 4); bits = rd16(f + 14);
            if (len > 16) fseek(fp, (long)(len - 16), SEEK_CUR);
        } else if (!memcmp(ch, "data", 4)) {
            size_t frames, i, bps = (size_t)bits / 8;
            unsigned char *raw;
            long here = ftell(fp), end;
            fseek(fp, 0, SEEK_END); end = ftell(fp); fseek(fp, here, SEEK_SET);
            if (len == 0 || len == 0xffffffffu || (long)len > end - here) len = (uint32_t)(end - here);
            if ((fmt != 1 && fmt != 3) || channels < 1 || (bits != 16 && bits != 8 && bits != 24 && bits != 32)) break;
            frames = len / (bps * (size_t)channels);
            raw = (unsigned char *)malloc(len);
            if (!raw || fread(raw, 1, len, fp) != len) { free(raw); break; }
            out = (float *)malloc(sizeof(float) * (frames ? frames : 1));
            for (i = 0; i < frames; i++) {
                const unsigned char *s = raw + i * bps * (size_t)channels;
                float v;
                if (fmt == 3 && bits == 32) { float f; memcpy(&f, s, 4); v = f; }
                else if (bits == 16) v = (float)(int16_t)rd16(s) / 32768.0f;
                else if (bits == 8) v = ((int)s[0] - 128) / 128.0f;
                else if (bits == 24) v = (float)((int32_t)((s[0] << 8) | (s[1] << 16) | ((uint32_t)s[2] << 24)) >> 8) / 8388608.0f;
                else v = (float)(int32_t)rd32(s) / 2147483648.0f;
                out[i] = v;
            }
            free(raw);
            *n = frames;
            break;
        } else {
            fseek(fp, (long)(len + (len & 1)), SEEK_CUR);
        }
    }
    fclose(fp);
    return out;
}

static void wr32(FILE *fp, uint32_t v) { unsigned char b[4] = { v & 255, (v >> 8) & 255, (v >> 16) & 255, (v >> 24) & 255 }; fwrite(b, 1, 4, fp); }
static void wr16(FILE *fp, uint16_t v) { unsigned char b[2] = { v & 255, (v >> 8) & 255 }; fwrite(b, 1, 2, fp); }

int wav_write(const char *path, const float *x, size_t n, int fs)
{
    FILE *fp = fopen(path, "wb");
    size_t i;
    if (!fp) return -1;
    fwrite("RIFF", 1, 4, fp); wr32(fp, (uint32_t)(36 + 2 * n)); fwrite("WAVE", 1, 4, fp);
    fwrite("fmt ", 1, 4, fp); wr32(fp, 16); wr16(fp, 1); wr16(fp, 1); wr32(fp, (uint32_t)fs);
    wr32(fp, (uint32_t)fs * 2); wr16(fp, 2); wr16(fp, 16);
    fwrite("data", 1, 4, fp); wr32(fp, (uint32_t)(2 * n));
    for (i = 0; i < n; i++) {
        float v = x[i] * 32767.0f;
        int s = (int)(v < 0 ? v - 0.5f : v + 0.5f);
        if (s > 32767) s = 32767;
        if (s < -32768) s = -32768;
        wr16(fp, (uint16_t)(int16_t)s);
    }
    fclose(fp);
    return 0;
}
