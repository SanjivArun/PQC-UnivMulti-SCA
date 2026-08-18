/*
 * Kyber Masked Decryption - Phases 2 & 3: Key-Sharing + Masked Dot Product
 *
 * Phase 2: arithmetic share splitting (coefficient domain)
 * Phase 3: masked dot product (NTT domain)
 *
 * Arithmetic sharing: s = s1 + s2 (mod q), x = x1 + x2 (mod q).
 * The secret key after polyvec_frombytes is in NTT domain.
 */

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "masked_poly.h"
#include "ntt.h"
#include "randombytes.h"

/* ======================== Phase 2: Key-Sharing ======================== */

static void random_mod_q(int16_t *out)
{
    uint8_t buf[2];
    do {
        buf[0] = 0;
        buf[1] = 0;
        randombytes(buf, 2);
        *out = ((buf[0] << 8) | buf[1]) & 0x0FFF;
    } while (*out >= KYBER_Q);
}

void masked_poly_split(poly *s1, poly *s2, const poly *s)
{
    unsigned int i;

    for (i = 0; i < KYBER_N; i++) {
        random_mod_q(&s1->coeffs[i]);
        s2->coeffs[i] = s->coeffs[i] - s1->coeffs[i];
        if (s2->coeffs[i] < 0) {
            s2->coeffs[i] += KYBER_Q;
        }
    }
}

/* ===================== Phase 3: Masked Dot Product ==================== */

/*
 * Masked dot product — share 1.
 * Computes out = -s1_ntt * b_ntt in NTT domain (pointwise basemul, accumulate, reduce).
 * Result stored in coefficient domain after implicit INTT by caller.
 *
 * The negation ensures x1 + x2 = -(s1+s2)*b = -(s*b).
 * After INTT: INTT(x1) + INTT(x2) = INTT(-s*b), so v - s*b = v + INTT(x1+x2).
 */
void masked_polyvec_basemul_acc_share1(poly *out, const polyvec *s1_ntt, const polyvec *b_ntt)
{
    int i;
    poly t;

    /* out = -s1_ntt[0] * b_ntt[0] */
    poly_basemul_montgomery(&t, &s1_ntt->vec[0], &b_ntt->vec[0]);
    for (i = 0; i < KYBER_N; i++) {
        out->coeffs[i] = -t.coeffs[i];
    }

    /* accumulate remaining: out += -(s1_ntt[i] * b_ntt[i]) */
    for (i = 1; i < KYBER_K; i++) {
        poly_basemul_montgomery(&t, &s1_ntt->vec[i], &b_ntt->vec[i]);
        poly_sub(out, out, &t);   /* out -= t => out += (-t) */
    }

    poly_reduce(out);
}

/*
 * Masked dot product — share 2.
 * Computes out = s2_ntt * b_ntt in NTT domain (pointwise basemul, accumulate, reduce).
 *
 * The v subtraction happens in Phase 5 after INTT: x2 ← intt(x2) - v.
 * This avoids doing NTT(v) and keeps the structure simpler.
 *
 * Combined with share1: x1_coeff + x2_coeff = -s*b in NTT domain.
 * After INTT: INTT(x1) + INTT(x2) = INTT(-s*b).
 * Then: v + INTT(x1+x2) = v - s*b ✓
 */
void masked_polyvec_basemul_acc_share2(poly *out, const polyvec *s2_ntt, const polyvec *b_ntt)
{
    int i;
    poly t;

    /* out = s2_ntt[0] * b_ntt[0] */
    poly_basemul_montgomery(out, &s2_ntt->vec[0], &b_ntt->vec[0]);

    /* accumulate: out += s2_ntt[i] * b_ntt[i] */
    for (i = 1; i < KYBER_K; i++) {
        poly_basemul_montgomery(&t, &s2_ntt->vec[i], &b_ntt->vec[i]);
        poly_add(out, out, &t);
    }

    poly_reduce(out);
}
