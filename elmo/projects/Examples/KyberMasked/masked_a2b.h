/*
 * Kyber Masked Decryption - Masked A2B compression (Phase 4)
 *
 * masked_poly_tomsg() computes the Boolean shares of Compress_q(x, 1)
 * from two ARITHMETIC shares of x without reconstructing x. It uses the
 * secure BDV21 single-lookup table-based A2B conversion
 * (ePrint 2021/067, Algorithm 8) together with the OSFP18 mod-q ->
 * mod-2^16 transform (no rejection sampling), exactly as in the
 * reference mkm4 implementation behind ePrint 2022/058.
 */

#ifndef MASKED_A2B_H
#define MASKED_A2B_H

#include <stdint.h>
#include "poly.h"

/*
 * Masked poly_tomsg.
 *
 * Converts two arithmetic shares x1, x2 (each mod-3329, with
 * x1[j] + x2[j] = x[j] mod 3329) to two Boolean shares m1, m2 such that
 * m1[j] XOR m2[j] = Compress_q(x[j], 1).
 *
 * The table generation (with fresh, per-call random masks) happens
 * internally — no external setup is needed.
 */
void masked_poly_tomsg(uint8_t m1[KYBER_INDCPA_MSGBYTES],
                       uint8_t m2[KYBER_INDCPA_MSGBYTES],
                       const poly *x1,
                       const poly *x2);

#endif
