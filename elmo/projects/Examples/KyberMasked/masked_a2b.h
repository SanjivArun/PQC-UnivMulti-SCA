/*
 * Kyber Masked Decryption - Masked A2B compression (Phase 4)
 *
 * This file provides the BDV21 "Debrāize (Fixed)" single-lookup
 * arithmetic-to-Boolean conversion for masked Kyber decryption,
 * plus the 4-step pipeline to convert arithmetic shares x1,x2 of
 * x ∈ Z_q into Boolean shares m1,m2 of Compress_q(x,1).
 *
 * See Sections 2.4 and 4 of the implementation plan for the full math.
 */

#ifndef MASKED_A2B_H
#define MASKED_A2B_H

#include <stdint.h>
#include "poly.h"

/*
 * State for BDV21 A2B tables and embedded randomness.
 *
 * Layout:
 *   tables[i][beta][a] — 16-bit entry for iteration i (k+1 = 9 bits used),
 *   beta ∈ {0, 1}(carry-in), a ∈ [0, 2^k).
 *   r[i] — the k-bit fresh random mask r_i for iteration i.
 */
typedef struct {
    uint16_t tables[2][2][256];  /* 2 iterations, 2 betas, 256 entries */
    uint8_t  r[2];            /* r_0, r_1 — 8 bits each */
} a2b_state;

/**
 * Generate BDV21 A2B tables (once per decryption call).
 *
 * Fills state->tables[] with random masking per the BDV21 spec (Algorithm 7).
 * Requires: state must point to at least a2b_state (zero-initialized or stack memory).
 */
void a2b_generate_tables(a2b_state *state);

/**
 * BDV21 fixed single-lookup A2B conversion for one 16-bit coefficient.
 *
 * Input:
 *   y1, y2 — 16-bit arithmetic shares with y1 + y2 = x mod 2^16
 *   state — pregenerated A2B tables and randomness
 *
 * Output:
 *   b[0], b[1] — two 8-bit Boolean share chunks.
 *   Combined: B = b[0] | (b[1] << 8) represents a Boolean sharing of (y1 XOR y2).
 *   After XORing with r = r_lo | (r_hi << 8): B XOR r = x mod 2^16.
 *
 * The MSB of (B XOR r) equals Compress_q(x, 1).
 */
void a2b_convert_16bit(uint16_t b[2], uint16_t y1, uint16_t y2,
                       const a2b_state *state);

/**
 * Masked poly_tomsg — 4-step pipeline on all 256 coefficients.
 *
 * Converts two arithmetic shares x1, x2 mod-3329 to two Boolean shares
 * m1, m2 such that m1[j] XOR m2[j] = Compress_q((x1[j]+x2[j]) mod 3329, 1).
 *
 * Pipeline per coefficient:
 *   1. y1 = x1[j] - 832   (subtract round(q/4))
 *   2. Modulus switch: reinterpret as 16-bit unsigned
 *   3. y1 = y1 - 1665     (subtract round(q/2))
 *   4. BDV21 A2B + MSB extraction → compressed bit
 *
 * m1 and m2 are output as byte arrays (32 bytes each).
 *
 * Table generation is called internally — no external a2b_state setup needed.
 */
void masked_poly_tomsg(uint8_t m1[KYBER_INDCPA_MSGBYTES],
                       uint8_t m2[KYBER_INDCPA_MSGBYTES],
                       const poly *x1,
                       const poly *x2);

#endif
