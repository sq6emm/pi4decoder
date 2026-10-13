/*
 * PI4 source coding, convolutional code, interleaver and sync vector.
 * Reference: OZ2M's PI4Ino encoder (BeaconPI4CW.ino) as ported in
 * mgmbeacon's PI4Encoder.cpp; constants checked against the worked example
 * on https://rudius.net/oz2m/ngnb/pi4.htm (OZ7IGY).
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "pi4.h"
#include <ctype.h>
#include <string.h>

const char pi4_alphabet[39] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /";

const uint8_t pi4_sync[PI4_NSYM] = {
    0,0,1,0,0,1,1,1,1,0,1,0,1,0,1,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,
    1,1,1,0,0,1,1,1,1,1,0,0,1,1,0,1,1,1,1,0,1,0,1,1,0,1,1,0,1,0,
    0,0,0,0,1,1,1,1,1,0,1,0,1,0,0,0,0,0,1,1,1,1,1,0,1,0,0,1,0,0,
    1,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,
    0,1,1,1,0,1,1,0,1,0,1,0,1,0,0,0,0,1,1,1,0,0,0,0,1,1,
};

/* Bit-reversal interleaver: codeword bit p goes to channel symbol
 * pi4_interleave[p]. Scan i = 0..255, r = bitrev8(i); whenever r < 146 the
 * next codeword bit is placed at r. */
const uint8_t pi4_interleave[PI4_NSYM] = {
    0,128,64,32,96,16,144,80,48,112,8,136,72,40,104,24,88,56,120,4,132,68,36,100,20,84,
    52,116,12,140,76,44,108,28,92,60,124,2,130,66,34,98,18,82,50,114,10,138,74,42,106,26,
    90,58,122,6,134,70,38,102,22,86,54,118,14,142,78,46,110,30,94,62,126,1,129,65,33,97,
    17,145,81,49,113,9,137,73,41,105,25,89,57,121,5,133,69,37,101,21,85,53,117,13,141,77,
    45,109,29,93,61,125,3,131,67,35,99,19,83,51,115,11,139,75,43,107,27,91,59,123,7,135,
    71,39,103,23,87,55,119,15,143,79,47,111,31,95,63,127,
};

static const int variant_k[4] = { 40, 80, 96, 120 };
static const char *const variant_name[4] = { "PI4", "PI4-80", "PI4-96", "PI4-120" };

int pi4_variant_k(pi4_variant v) { return variant_k[v & 3]; }
double pi4_spacing_hz(pi4_variant v) { return variant_k[v & 3] * 12000.0 / 2048.0; }
const char *pi4_variant_name(pi4_variant v) { return variant_name[v & 3]; }

int pi4_variant_parse(const char *s)
{
    int i;
    for (i = 0; i < 4; i++) {
        const char *a = variant_name[i], *b = s;
        while (*a && *b && toupper((unsigned char)*a) == toupper((unsigned char)*b)) a++, b++;
        if (!*a && !*b) return i;
    }
    if (!strcmp(s, "80") || !strcmp(s, "96") || !strcmp(s, "120"))
        return s[0] == '8' ? 1 : s[0] == '9' ? 2 : 3;
    return -1;
}

int pi4_pack(const char *msg, uint64_t *n)
{
    uint64_t v = 0;
    int i;
    for (i = 0; i < PI4_NCHAR; i++) {
        int c = *msg ? toupper((unsigned char)*msg++) : ' ';
        const char *p = c ? strchr(pi4_alphabet, c) : NULL;
        if (!p) return -1;
        v = v * 38 + (uint64_t)(p - pi4_alphabet);
    }
    if (*msg) return -1;
    *n = v;
    return 0;
}

void pi4_unpack(uint64_t n, char out[PI4_NCHAR + 1])
{
    int i;
    out[PI4_NCHAR] = 0;
    for (i = PI4_NCHAR - 1; i >= 0; i--) {
        out[i] = pi4_alphabet[n % 38];
        n /= 38;
    }
    for (i = PI4_NCHAR - 1; i >= 0 && out[i] == ' '; i--) out[i] = 0;
}

static int parity32(uint32_t x)
{
    x ^= x >> 16; x ^= x >> 8; x ^= x >> 4; x ^= x >> 2; x ^= x >> 1;
    return (int)(x & 1);
}

void pi4_coded_bits(uint64_t n, uint8_t coded[PI4_NSYM])
{
    uint32_t reg = 0;
    int j, idx = 0;
    for (j = 0; j < PI4_NBITS; j++) {
        reg = (reg << 1) | (uint32_t)((n >> 41) & 1);
        n <<= 1;
        coded[idx++] = (uint8_t)parity32(reg & 0xF2D05351u);
        coded[idx++] = (uint8_t)parity32(reg & 0xE4613C47u);
    }
}

void pi4_encode_n(uint64_t n, uint8_t tones[PI4_NSYM])
{
    uint8_t coded[PI4_NSYM];
    int i;
    pi4_coded_bits(n, coded);
    for (i = 0; i < PI4_NSYM; i++) tones[pi4_interleave[i]] = coded[i];
    for (i = 0; i < PI4_NSYM; i++) tones[i] = (uint8_t)(pi4_sync[i] + 2 * tones[i]);
}

int pi4_encode(const char *msg, uint8_t tones[PI4_NSYM])
{
    uint64_t n;
    if (pi4_pack(msg, &n)) return -1;
    pi4_encode_n(n, tones);
    return 0;
}

static int is_alpha(char c) { return c >= 'A' && c <= 'Z'; }
static int is_digit(char c) { return c >= '0' && c <= '9'; }

/* Amateur callsign shape: prefix ([A-Z]{1,2} or [0-9][A-Z]{1,2}), 1-2
 * digits, suffix [A-Z]{1,3}. Covers OZ7IGY, ED5YAE, IW9GDC, G3ZJO, SR6LEG,
 * 4X1AB, 9A1AA, 3DA0AB; rejects random strings such as 8DN1FWIP. */
static int callsign_like(const char *s, size_t n)
{
    size_t i = 0, k;
    if (n < 3 || n > 7) return 0;
    if (is_digit(s[0])) {
        i = 1;
        for (k = 0; i < n && is_alpha(s[i]) && k < 2; i++, k++) ;
        if (k == 0) return 0;
    } else {
        for (k = 0; i < n && is_alpha(s[i]) && k < 2; i++, k++) ;
        if (k == 0) return 0;
    }
    for (k = 0; i < n && is_digit(s[i]) && k < 2; i++, k++) ;
    if (k == 0) return 0;
    for (k = 0; i < n && is_alpha(s[i]) && k < 3; i++, k++) ;
    if (k == 0) return 0;
    return i == n;
}

/* A PI4 message is normally the beacon callsign, optionally with a "/B"
 * style suffix or a "OZ/" style prefix; the protocol also allows status
 * strings starting with '/'. Anything else is "implausible" (still
 * reported, but only through a stricter gate). */
int pi4_text_plausible(const char *text)
{
    size_t len = strlen(text);
    const char *slash = strchr(text, '/');
    if (len < 2 || len > PI4_NCHAR) return 0;
    if (text[0] == '/') return !strchr(text + 1, ' ') && !strchr(text + 1, '/');
    if (strchr(text, ' ')) return 0;
    if (!slash) return callsign_like(text, len);
    if (strchr(slash + 1, '/')) return 0;
    if (callsign_like(text, (size_t)(slash - text))) {       /* CALL/B, CALL/P */
        size_t suf = len - (size_t)(slash - text) - 1, i;
        if (suf < 1 || suf > 2) return 0;
        for (i = 0; i < suf; i++) if (!is_alpha(slash[1 + i]) && !is_digit(slash[1 + i])) return 0;
        return 1;
    }
    if (slash - text >= 1 && slash - text <= 3) {            /* PREFIX/CALL */
        size_t i;
        for (i = 0; i < (size_t)(slash - text); i++) if (!is_alpha(text[i]) && !is_digit(text[i])) return 0;
        return callsign_like(slash + 1, len - (size_t)(slash - text) - 1);
    }
    return 0;
}
