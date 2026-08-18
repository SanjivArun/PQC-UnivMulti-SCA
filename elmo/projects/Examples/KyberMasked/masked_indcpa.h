/*
 * Kyber Masked Decryption - Phase 5: masked_indcpa_dec API
 *
 * Implements masked_indcpa_dec() — a full masked version of indcpa_dec()
 * that uses BDV21 A2B for the nonlinear compression step.
 *
 * Same signature as indcpa_dec() for drop-in compatibility.
 * Use #define KYBER_MASKED to replace indcpa_dec with masked_indcpa_dec.
 */

#ifndef MASKED_INDCPA_H
#define MASKED_INDCPA_H

#include <stdint.h>
#include "params.h"

/**
 * Masked Kyber IND-CPA decryption.
 *
 * Input:
 *   c — ciphertext (pack(b) + pack(v))
 *   sk — secret key (pack(s) where s is in NTT domain)
 *
 * Output:
 *   m — 32-byte message (identical output format to indcpa_dec)
 *
 * Internally splits secret key into arithmetic shares, computes masked
 * dot product, then applies BDV21 A2B for masked compression.
 * The secret key s is never reconstructed in full during decryption.
 */
void masked_indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                      const uint8_t c[KYBER_INDCPA_BYTES],
                      const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

#endif
