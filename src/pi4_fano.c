/*
 * Sequential (Fano) decoder for the rate-1/2 K=32 convolutional code
 * (Layland-Lushbaugh polynomials 0xF2D05351 / 0xE4613C47, the code shared
 * by WSPR, JT4 and PI4). Written from the textbook algorithm; integer
 * metrics, no allocation, 73 nodes of state on the stack.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "pi4.h"

#define POLY1 0xF2D05351u
#define POLY2 0xE4613C47u

static uint32_t bitrev32(uint32_t x)
{
    uint32_t r = 0; int i;
    for (i = 0; i < 32; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}

struct node {
    uint32_t state;   /* encoder register after this node's decision        */
    long     gamma;   /* path metric up to this node                        */
    int      tm[2];   /* branch metrics, best first                         */
    uint8_t  bit[2];  /* input bit of each branch                           */
    uint8_t  i;       /* branch currently being tried                       */
};

static int par(uint32_t x)
{
    x ^= x >> 16; x ^= x >> 8; x ^= x >> 4; x ^= x >> 2; x ^= x >> 1;
    return (int)(x & 1);
}

/* Both polynomials have bit 0 set, so flipping the input bit flips both
 * coded bits: sym(bit=1) == sym(bit=0) ^ 3. */
static void expand(struct node *np, const int16_t *metric, int j, uint32_t prev, int rev)
{
    uint32_t s0 = prev << 1, p1 = rev ? bitrev32(POLY1) : POLY1, p2 = rev ? bitrev32(POLY2) : POLY2;
    int pair = rev ? PI4_NBITS - 1 - j : j;
    int sym0 = (par(s0 & p1) << 1) | par(s0 & p2);
    int m0 = metric[2 * (2 * pair) + (sym0 >> 1)] + metric[2 * (2 * pair + 1) + (sym0 & 1)];
    np->i = 0;
    if (j >= PI4_NINFO) {                  /* tail: zero branch only */
        np->tm[0] = m0; np->bit[0] = 0;
        np->tm[1] = -1000000; np->bit[1] = 0;
        return;
    }
    {
        int sym1 = sym0 ^ 3;
        int m1 = metric[2 * (2 * pair) + (sym1 >> 1)] + metric[2 * (2 * pair + 1) + (sym1 & 1)];
        if (m0 >= m1) { np->tm[0] = m0; np->bit[0] = 0; np->tm[1] = m1; np->bit[1] = 1; }
        else          { np->tm[0] = m1; np->bit[0] = 1; np->tm[1] = m0; np->bit[1] = 0; }
    }
}

int pi4_fano(const int16_t *metric, int delta, long maxcycles, int reverse,
             uint8_t bits[PI4_NBITS], long *metric_out, long *cycles_out)
{
    struct node nodes[PI4_NBITS + 1];
    struct node *np = nodes, *last = nodes + PI4_NBITS;
    long t = 0, cycles = 0;
    int j;

    np->state = 0; np->gamma = 0;
    expand(np, metric, 0, 0, reverse);
    for (;;) {
        long ngamma = np->gamma + np->tm[np->i];
        j = (int)(np - nodes);
        if (ngamma >= t) {
            if (np->gamma < t + delta)         /* first visit: tighten */
                while (ngamma >= t + delta) t += delta;
            np[1].gamma = ngamma;
            np[1].state = (np->state << 1) | np->bit[np->i];
            np++;
            if (np == last) break;             /* all 73 bits decided */
            expand(np, metric, j + 1, np->state, reverse);
        } else {
            for (;;) {                          /* back up / lateral */
                if (np == nodes || np[-1].gamma < t) {
                    t -= delta;                 /* cannot back up: lower T */
                    np->i = 0;
                    break;
                }
                np--;
                if ((int)(np - nodes) < PI4_NINFO && np->i == 0) {
                    np->i = 1;                  /* try the other branch */
                    break;
                }
            }
            if (++cycles > maxcycles) {
                if (cycles_out) *cycles_out = cycles;
                if (metric_out) *metric_out = np->gamma;
                return 0;
            }
        }
    }
    for (j = 0; j < PI4_NBITS; j++) bits[j] = (uint8_t)(nodes[j + 1].state & 1);
    if (reverse) {                         /* info bits came out last-first */
        for (j = 0; j < PI4_NINFO / 2; j++) { uint8_t t = bits[j]; bits[j] = bits[PI4_NINFO - 1 - j]; bits[PI4_NINFO - 1 - j] = t; }
    }
    if (metric_out) *metric_out = last->gamma;
    if (cycles_out) *cycles_out = cycles;
    return 1;
}
