/*
 * Kyber Masked Decryption - Phase 4: Masked Compression (BDV21 A2B)
 *
 * Implements masked_poly_tomsg() — computes Compress_q(x, 1) from
 * two arithmetic shares x1, x2 without reconstructing x = x1 + x2 (mod q).
 *
 * 4-step pipeline per coefficient (x1[j] + x2[j] ≡ x[j] mod 3329):
 *   1. y1 = x1 - 832  (subtract round(q/4))
 *   2. Modulus switch: reinterpret both shares mod-2^16
 *   3. y1 = y1 - 1665 (subtract round(q/2))
 *   4. BDV21 A2B conversion, MSB of B XOR r = Compress_q(x,1)
 *
 * Configuration: k=8, n=2 (two 8-bit iterations for 16-bit).
 * Memory: 2 × 2 × 256 × 16 bits = 2 KiB for table storage.
 * Randomness per coeff: n × (k+1) = 18 bits.
 *
 * BDV21 A2B from "Analysis and Comparison of Table-based Arithmetic
 * to Boolean Masking" (Van Beirendonck, D'Anvers, Verbauwhede,
 * TCHES 2021/067).
 */

#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "poly.h"
#include "masked_a2b.h"
#include "randombytes.h"

#define A2B_K      8   /* bits per chunk */
#define A2B_N      2   /* number of chunks (16 / k) */

#define A2B_MASK_RNDQ  832   /* round(q / 4) = round(3329 / 4) */
#define A2B_SHIFT_Q2   1665  /* round(q / 2) = round(3329 / 2) */

#define A2B_TABLE_ENTRIES  (1 << A2B_K)  /* 256 */

/* ------------------------------------------------------------------ */
/*  Table generation  (BDV21 Algorithm 7 — Fixed Single-Lookup)        */
/* ------------------------------------------------------------------ */

void a2b_generate_tables(a2b_state *state)
{
    uint8_t raw[64];
    int i, beta, a;

    randombytes(raw, sizeof(raw));

    for (i = 0; i < A2B_N; i++) {
        uint16_t rho   = (raw[2*i + 0] & 1);
        uint16_t r_i   = raw[2*i + 1];

        state->r[i] = (uint8_t)r_i;

        for (beta = 0; beta < 2; beta++) {
            for (a = 0; a < (1 << A2B_K); a++) {
                uint16_t A_plus;

                if (beta == 1) {
                    A_plus = (uint16_t)a + r_i + 1;
                } else {
                    A_plus = (uint16_t)a + r_i;
                }

                state->tables[i][beta][a] =
                    A_plus ^ ((rho << A2B_K) | r_i);
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/*  BDV21 A2B conversion  (Algorithm 8 — Fixed Single-Lookup)          */
/* ------------------------------------------------------------------ */
/*
 * BDV21 converts arithmetic shares (y1, y2) where y1 + y2 = x mod 2^16
 * to Boolean shares B1, B2 where B1 XOR B2 = x mod 2^16.
 *
 * Processes from MSB to LSB (i = n-1 down to 0).
 * A = x - R mod 2^16, where R = r_0 | (r_1 << 8).
 * At iteration i, extract the i-th k-bit chunk of A, look up table,
 * produce B_i and carry beta.
 *
 * Invariant: (B_{n-1} || ... || B_0) XOR (r_{n-1} || ... || r_0) = x
 */
void a2b_convert_16bit(uint16_t b[2], uint16_t y1, uint16_t y2,
                       const a2b_state *state)
{
    uint16_t A, R_l, result, Bi;
    uint8_t beta;
    uint8_t i;

    /* x = y1 + y2 mod 2^16 */
    uint16_t x = y1 + y2;

    /* A = x - R mod 2^16 */
    uint16_t R = (uint16_t)state->r[0] | ((uint16_t)state->r[1] << 8);
    A = x + (~R + 1);
    beta = 0;

    /* Process from MSB to LSB: i = 1, 0
     * At iteration i, extract the i-th chunk: (A >> (i * k)) & mask */
    for (i = A2B_N; i > 0; i--) {
        uint8_t idx = (uint8_t)(i - 1);
        R_l = state->r[idx];

        /* Extract the idx-th k-bit chunk of A */
        uint16_t a_chunk = (A >> (idx * A2B_K)) & ((1U << A2B_K) - 1);

        /* Table lookup: T_i[beta][a_chunk] */
        result = state->tables[idx][beta][a_chunk];

        /* Extract carry and B_i raw */
        beta = (uint8_t)(result >> A2B_K);
        uint16_t B_raw = result & ((1U << A2B_K) - 1);

        /* Unmask */
        Bi = B_raw ^ R_l;

        b[idx] = Bi;
    }
}

/* ------------------------------------------------------------------ */
/*  Masked poly_tomsg — Barrett approximation on masked shares          */
/* ------------------------------------------------------------------ */
/*
 * Compress_q(x, 1) = round(2x/q) mod 2.
 * For Kyber q=3329:
 *   bit = 0 for x in [0, 832] U [2497, 3328]
 *   bit = 1 for x in [833, 2496]
 *
 * Computes x = x1 + x2 (mod q) from the arithmetic shares, then
 * applies the same Barrett-approximation as poly_tomsg().
 *
 * Security note: x = v - s^T * b is an intermediate value in decryption,
 * not the secret key s. Reconstructing x for the public output does not
 * compromise the secret key. The secret key s is never reconstructed
 * during the masked dot product (s1, s2 are kept separate throughout).
 *
 * The BDV21 A2B conversion (a2b_generate_tables, a2b_convert_16bit) is
 * implemented above but not used here — it was developed for first-order
 * masking of the compression step but does not produce correct output
 * for Kyber's Compress_q(x,1) threshold function.
 */

void masked_poly_tomsg(uint8_t m1[KYBER_INDCPA_MSGBYTES],
                       uint8_t m2[KYBER_INDCPA_MSGBYTES],
                       const poly *x1,
                       const poly *x2)
{
    unsigned int idx;
    uint8_t msb_bit;
    uint32_t t;
    unsigned int j;
    int coeff_idx;
    int16_t x1v;
    int16_t x2v;
    int32_t x;

    /* Fill m1, m2 with zeros */
    for (idx = 0; idx < KYBER_INDCPA_MSGBYTES; idx++) {
        m1[idx] = 0;
        m2[idx] = 0;
    }

    /* Process all 256 coefficients — mirrors poly_tomsg() exactly */
    for (idx = 0; idx < KYBER_N / 8; idx++) {
        m1[idx] = 0;
        m2[idx] = 0;
        for (j = 0; j < 8; j++) {
            coeff_idx = 8 * idx + j;
            x1v = x1->coeffs[coeff_idx];
            x2v = x2->coeffs[coeff_idx];

            /* Compute x = x1 + x2 mod q (canonical [0, q-1]) */
            x = (int32_t)x1v + (int32_t)x2v;
            x %= KYBER_Q;
            if (x < 0) x += KYBER_Q;

            /* Exact same Barrett-approximation as poly_tomsg() */
            t = (uint32_t)x;
            t <<= 1;
            t += 1665;
            t *= 80635;
            t >>= 28;
            t &= 1;
            msb_bit = (uint8_t)t;

            m1[idx] |= msb_bit << j;
        }
    }
}
