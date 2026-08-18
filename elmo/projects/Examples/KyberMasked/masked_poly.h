/*
 * Kyber Masked Decryption - Phase 2: Key-Sharing Infrastructure (Public Header)
 *
 * Arithmetic share splitting for Kyber polynomials and polyvecs.
 */

#ifndef MASKED_POLY_H
#define MASKED_POLY_H

#include <stdint.h>
#include "poly.h"
#include "polyvec.h"

/*
 * Split polynomial s into two arithmetic shares s1, s2.
 * s = s1 + s2 (mod q).
 * s1 is uniformly random in [0, q-1].
 * s2 = s - s1 (mod q), also in [0, q-1].
 * All inputs/outputs in coefficient domain (not NTT domain).
 * Each coefficient must be in [0, q-1] on input.
 */
void masked_poly_split(poly *s1, poly *s2, const poly *s);

/*
 * Masked dot product — share 1.
 * Computes out = -s1_ntt * b_ntt via NTT-domain basemul.
 * Result is in coefficient domain after caller applies INTT.
 * The negation ensures x1 + x2 = -(s1+s2)*b = -s*b.
 * s1_ntt and b_ntt must both be in NTT domain.
 */
void masked_polyvec_basemul_acc_share1(poly *out, const polyvec *s1_ntt, const polyvec *b_ntt);

/*
 * Masked dot product — share 2.
 * Computes out = s2_ntt * b_ntt via NTT-domain basemul.
 * Result is in coefficient domain after caller applies INTT.
 * The v subtraction happens in Phase 5 after INTT: x2 ← intt(x2) - v.
 * s2_ntt and b_ntt must both be in NTT domain.
 */
void masked_polyvec_basemul_acc_share2(poly *out, const polyvec *s2_ntt, const polyvec *b_ntt);

#endif
