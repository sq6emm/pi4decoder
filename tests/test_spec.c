/* Spec round-trip tests against the OZ2M worked example (OZ7IGY) and the
 * reference encoder's construction. SPDX-License-Identifier: GPL-3.0-or-later */
#include "pi4.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(c, ...) do { if (!(c)) { fails++; printf("FAIL %s:%d ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

/* OZ7IGY symbols from https://rudius.net/oz2m/ngnb/pi4.htm (also the output
 * of OZ2M's PI4Ino reference encoder). */
static const uint8_t oz7igy[PI4_NSYM] = {
    2,0,1,0,0,3,3,3,3,2,3,2,1,2,1,2,0,3,2,2,0,3,2,2,0,1,1,0,0,1,
    3,1,3,0,2,1,1,3,3,1,2,0,1,3,2,1,3,3,3,2,1,2,3,1,2,1,1,0,3,2,
    0,2,0,0,1,3,3,1,3,2,3,2,3,0,2,0,0,2,1,3,3,3,1,2,3,0,0,3,0,2,
    3,2,1,0,2,0,2,1,0,0,1,1,0,2,0,2,2,3,3,2,2,2,2,3,1,0,0,1,3,3,
    0,1,3,1,2,1,3,0,3,0,3,0,1,2,2,0,2,3,1,3,2,0,0,2,1,1,
};

/* reference encoder transcribed from PI4Encoder.cpp (OZ2M PI4Ino port) */
static int parity32(uint32_t v) { int e = 0, b; for (b = 0; b < 32; b++) if ((v >> b) & 1) e = 1 - e; return e; }
static void ref_encode(const char *message, uint8_t tones[PI4_NSYM])
{
    const char *chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /";
    uint64_t src = 0; uint32_t n = 0; uint8_t conv[PI4_NSYM], il[PI4_NSYM]; int i, j, idx = 0, p = 0;
    char padded[9]; memset(padded, ' ', 8); padded[8] = 0; memcpy(padded, message, strlen(message));
    for (i = 0; i < 8; i++) {
        char c = padded[i];
        src = src * 38 + (uint64_t)(strchr(chars, c) - chars);
    }
    for (j = 0; j < PI4_NSYM / 2; j++) {
        n <<= 1;
        if (src & 0x20000000000ULL) n |= 1;
        src <<= 1;
        conv[idx++] = (uint8_t)parity32(n & 0xF2D05351);
        conv[idx++] = (uint8_t)parity32(n & 0xE4613C47);
    }
    memset(il, 0, sizeof il);
    for (i = 0; i <= 255; i++) {
        int r = 0, b;
        for (b = 0; b <= 7; b++) if ((i >> b) & 1) r |= 1 << (7 - b);
        if (p < PI4_NSYM && r < PI4_NSYM) il[r] = conv[p++];
    }
    for (i = 0; i < PI4_NSYM; i++) tones[i] = (uint8_t)(pi4_sync[i] + 2 * il[i]);
}

int main(void)
{
    uint8_t t[PI4_NSYM], r[PI4_NSYM];
    uint64_t n;
    char txt[9];
    int i;
    const char *msgs[] = { "OZ7IGY", "SR6LEG", "ED5YAE/B", "IW9GDC/B", "G3ZJO", "PE1ITR/B", "/ABCDEFG", "A", "00000000", "ZZZZZZZZ", "SQ6EMM" };

    /* worked example */
    CHECK(pi4_pack("OZ7IGY", &n) == 0 && n == 2851949862724ULL, "pack OZ7IGY = %llu", (unsigned long long)n);
    CHECK(pi4_encode("OZ7IGY", t) == 0 && !memcmp(t, oz7igy, PI4_NSYM), "OZ7IGY symbols");
    pi4_unpack(n, txt); CHECK(!strcmp(txt, "OZ7IGY"), "unpack -> '%s'", txt);

    /* interleaver table equals the spec's construction */
    {   int p = 0; uint8_t il[PI4_NSYM];
        for (i = 0; i <= 255; i++) {
            int rr = 0, b;
            for (b = 0; b <= 7; b++) if ((i >> b) & 1) rr |= 1 << (7 - b);
            if (p < PI4_NSYM && rr < PI4_NSYM) il[p++] = (uint8_t)rr;
        }
        CHECK(!memcmp(il, pi4_interleave, PI4_NSYM), "interleave table");
    }

    /* round trips against the reference encoder */
    for (i = 0; i < (int)(sizeof msgs / sizeof msgs[0]); i++) {
        CHECK(pi4_encode(msgs[i], t) == 0, "encode %s", msgs[i]);
        ref_encode(msgs[i], r);
        CHECK(!memcmp(t, r, PI4_NSYM), "reference encoder agrees for %s", msgs[i]);
        pi4_pack(msgs[i], &n); pi4_unpack(n, txt);
        CHECK(!strcmp(txt, msgs[i]), "round trip %s -> %s", msgs[i], txt);
    }
    CHECK(pi4_pack("TOOLONG123", &n) != 0, "reject long");
    CHECK(pi4_pack("OZ7IGY!", &n) != 0, "reject bad char");
    CHECK(pi4_pack("oz7igy", &n) == 0, "lower case accepted");

    /* sync bit is the low bit of every symbol */
    for (i = 0; i < PI4_NSYM; i++) CHECK((oz7igy[i] & 1) == pi4_sync[i], "sync bit %d", i);

    /* Fano decodes noiseless metrics of every message */
    for (i = 0; i < (int)(sizeof msgs / sizeof msgs[0]); i++) {
        int16_t met[PI4_NSYM][2]; uint8_t bits[PI4_NBITS]; int k; uint64_t m = 0; long metric, cyc;
        pi4_encode(msgs[i], t);
        for (k = 0; k < PI4_NSYM; k++) {
            int bit = t[pi4_interleave[k]] >> 1;
            met[k][bit] = 16; met[k][!bit] = -200;
        }
        CHECK(pi4_fano(&met[0][0], 64, 10000, 0, bits, &metric, &cyc) == 1 && pi4_fano(&met[0][0], 64, 10000, 1, bits, &metric, &cyc) == 1, "fano converges %s", msgs[i]);
        for (k = 0; k < PI4_NINFO; k++) m = (m << 1) | bits[k];
        pi4_pack(msgs[i], &n);
        CHECK(m == n, "fano recovers %s (got %llu want %llu)", msgs[i], (unsigned long long)m, (unsigned long long)n);
        for (k = PI4_NINFO; k < PI4_NBITS; k++) CHECK(bits[k] == 0, "tail zero");
    }

    /* text plausibility */
    {   const char *good[] = { "OZ7IGY", "ED5YAE/B", "IW9GDC/B", "G3ZJO", "SR6LEG", "4X1AB", "9A1AA", "3DA0AB", "/STATUS", "OZ/DL1AB", "PE1ITR/P", "SK6YH" };
        const char *bad[] = { "8DN1FWIP", "00000000", "AAAAAAAA", "VUOJEZIQ", "A", "OZ7IGY X", "12345", "ABCDEFGH", "OZ7IGY/", "/" };
        unsigned k;
        for (k = 0; k < sizeof good / sizeof good[0]; k++) CHECK(pi4_text_plausible(good[k]), "plausible %s", good[k]);
        for (k = 0; k < sizeof bad / sizeof bad[0]; k++) CHECK(!pi4_text_plausible(bad[k]), "implausible %s", bad[k]);
    }

    CHECK(pi4_spacing_hz(PI4_V_PI4) == 234.375 && pi4_spacing_hz(PI4_V_80) == 468.75 &&
          pi4_spacing_hz(PI4_V_96) == 562.5 && pi4_spacing_hz(PI4_V_120) == 703.125, "spacings");
    CHECK(pi4_variant_parse("pi4-96") == 2 && pi4_variant_parse("PI4") == 0 && pi4_variant_parse("x") < 0, "variant names");

    printf(fails ? "test_spec: %d FAILURES\n" : "test_spec: all passed\n", fails);
    return fails != 0;
}
