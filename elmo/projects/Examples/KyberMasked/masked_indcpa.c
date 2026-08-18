/*
 * Kyber Masked Decryption - Phase 5: Masked indcpa_dec
 *
 * Mirrors indcpa_dec() from indcpa.c but uses arithmetic shares
 * for all sensitive intermediates. The secret key s is split into
 * s1, s2 where s = s1 + s2 (mod q). The dot product is computed
 * as two masked shares x1 = -(s1^T * b) and x2 = v - (s2^T * b)
 * so that x1 + x2 = v - s^T * b without ever reconstructing s.
 *
 * The nonlinear compression step (poly_tomsg) is replaced with
 * masked_poly_tomsg() which uses BDV21 A2B to produce Boolean
 * shares m1, m2 of Compress_q(x, 1) from arithmetic shares x1, x2.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "masked_poly.h"
#include "masked_a2b.h"
#include "masked_indcpa.h"

void masked_indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                       const uint8_t c[KYBER_INDCPA_BYTES],
                       const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
    polyvec b, skpv;
    polyvec s1_ntt, s2_ntt;
    poly v, share1, share2;
    poly x1, x2;
    uint8_t msg1[KYBER_INDCPA_MSGBYTES], msg2[KYBER_INDCPA_MSGBYTES];
    unsigned int j;

    /* Step 1: Unpack ciphertext — same as indcpa.c:44-48 */
    polyvec_decompress(&b, c);
    poly_decompress(&v, c + KYBER_POLYVECCOMPRESSEDBYTES);

    /* Step 2: Unpack secret key — same as indcpa.c:33-36 */
    /* skpv is in NTT domain (Montgomery representation) */
    polyvec_frombytes(&skpv, sk);

    /* Step 3: Split secret key into arithmetic shares */
    /* s = s1 + s2 (mod q), both in NTT domain */
    for (j = 0; j < KYBER_K; j++) {
        masked_poly_split(&s1_ntt.vec[j], &s2_ntt.vec[j], &skpv.vec[j]);
    }

    /* Step 4: NTT ciphertext — same as indcpa.c:190 */
    polyvec_ntt(&b);

    /* Step 5: Masked dot product — replaces indcpa.c:191 */
    /* share1 = -(s1_ntt * b_ntt), share2 = (s2_ntt * b_ntt) */
    /* share1 + share2 = -(s1 + s2)_ntt * b_ntt = -s_ntt * b_ntt */
    masked_polyvec_basemul_acc_share1(&share1, &s1_ntt, &b);
    masked_polyvec_basemul_acc_share2(&share2, &s2_ntt, &b);

    /* Step 6: INTT to coefficient domain — same as indcpa.c:192 */
    poly_invntt_tomont(&share1);
    poly_invntt_tomont(&share2);

    /* Step 7: Build x1, x2 shares of (v - s^T * b) */
    /* x1 = share1 = -(s1.b), x2 = v - share2 = v - (s2.b) */
    /* x1 + x2 = v - (s1 + s2).b = v - s.b  (never reconstructed) */
    for (j = 0; j < KYBER_N; j++) {
        x1.coeffs[j] = share1.coeffs[j];
    }
    poly_sub(&x2, &v, &share2);

    /* Step 8: Modular reduction — same as indcpa.c:195 */
    poly_reduce(&x1);
    poly_reduce(&x2);

    /* Step 9: Masked message decode — replaces indcpa.c:197 */
    /* poly_tomsg(m, &mp) becomes masked_poly_tomsg(msg1, msg2, &x1, &x2) */
    /* Produces Boolean shares m1, m2 where m1[j] XOR m2[j] = Compress_q(x1[j]+x2[j], 1) */
    masked_poly_tomsg(msg1, msg2, &x1, &x2);

    /* Step 10: Reconstruct message — boolean shares XOR */
    for (j = 0; j < KYBER_INDCPA_MSGBYTES; j++) {
        m[j] = msg1[j] ^ msg2[j];
    }
}
