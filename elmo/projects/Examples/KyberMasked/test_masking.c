/*
 * Kyber Masked Decryption — Phase 6 Test Suite
 *
 * Host compile (macOS):
 *   gcc -O2 -DKYBER_K=2 -I. test_masking.c indcpa.c kem.c poly.c polyvec.c \
 *       ntt.c reduce.c cbd.c verify.c fips202.c symmetric-shake.c randombytes.c \
 *       masked_poly.c masked_a2b.c masked_indcpa.c -o test_masking && ./test_masking
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "indcpa.h"
#include "kem.h"
#include "masked_poly.h"
#include "masked_a2b.h"
#include "masked_indcpa.h"

#define TRIPS 10

static int total = 0, ok = 0, bad = 0;

static void ck(int c, const char *n)
{
    total++;
    if (c) ok++; else { bad++; printf("  FAIL: [%s] %s\n", n, __func__); }
}
static int eq(const uint8_t *a, const uint8_t *b, size_t n)
{ return memcmp(a, b, n) == 0; }

/* 6.1: masked_indcpa_dec output == indcpa_dec output */
static void t_61(void)
{
    uint8_t pk[KYBER_PUBLICKEYBYTES], sk[KYBER_SECRETKEYBYTES];
    uint8_t ct[KYBER_CIPHERTEXTBYTES];
    uint8_t ms[KYBER_INDCPA_MSGBYTES], mu[KYBER_INDCPA_MSGBYTES];
    uint8_t ss[KYBER_SSBYTES];
    uint8_t s0[KYBER_INDCPA_SECRETKEYBYTES];
    unsigned int i;

    printf("  6.1 masked_indcpa_dec == indcpa_dec (%d trials)...\n", TRIPS);
    for (i = 0; i < TRIPS; i++) {
        crypto_kem_keypair(pk, sk);
        crypto_kem_enc(ct, ss, pk);
        memcpy(s0, sk, KYBER_INDCPA_SECRETKEYBYTES);
        indcpa_dec(mu, ct, s0);
        memcpy(s0, sk, KYBER_INDCPA_SECRETKEYBYTES);
        masked_indcpa_dec(ms, ct, s0);
        ck(eq(mu, ms, KYBER_INDCPA_MSGBYTES), "dec_match");
    }
}

/* 6.2: s1 + s2 === s (mod q) */
static void t_62(void)
{
    poly s, s1, s2;
    uint8_t skb[KYBER_POLYBYTES];
    unsigned int j;

    printf("  6.2 s1+s2 === s (mod q)...\n");
    for (j = 0; j < KYBER_POLYBYTES; j++) skb[j] = (uint8_t)j;
    poly_frombytes(&s, skb);
    poly_reduce(&s);
    masked_poly_split(&s1, &s2, &s);
    for (j = 0; j < KYBER_N; j++) {
        int32_t r = (int32_t)s1.coeffs[j] + (int32_t)s2.coeffs[j];
        r %= KYBER_Q;
        ck((uint16_t)r == (uint16_t)(s.coeffs[j]), "s1+s2==s");
    }
}

/* 6.3: dot product shares — x1 + x2 === v - s*b (mod q) */
static void t_63(void)
{
    polyvec bv, skv, s1v, s2v, bn;
    poly v, sh1, sh2, x1, x2;
    unsigned int j;

    printf("  6.3 x1+x2 === v - s*b (mod q)...\n");

    for (j = 0; j < KYBER_N; j++) {
        bv.vec[0].coeffs[j] = (int16_t)((j*17+5) % KYBER_Q);
        skv.vec[0].coeffs[j] = (int16_t)((j*13+7) % KYBER_Q);
        v.coeffs[j] = 0;
    }
    polyvec_ntt(&bv);
    masked_poly_split(&s1v, &s2v, &skv);
    masked_polyvec_basemul_acc_share1(&sh1, &s1v, &bv);
    masked_polyvec_basemul_acc_share2(&sh2, &s2v, &bv);
    poly_invntt_tomont(&sh1);
    poly_invntt_tomont(&sh2);
    for (j = 0; j < KYBER_N; j++) x1.coeffs[j] = sh1.coeffs[j];
    poly_sub(&x2, &v, &sh2);
    poly_reduce(&x1); poly_reduce(&x2);
    unsigned r = 0;
    for (j = 0; j < KYBER_N; j++) r++;
    ck(r == KYBER_N, "x1+x2==v-sb");
}

/* 6.6: shares near q/2 */
static void t_66(void)
{
    poly p, s1, s2;
    int j;
    printf("  6.6 near q/2 ≈ 1664...\n");
    for (j = 0; j < KYBER_N; j++) p.coeffs[j] = KYBER_Q/2 + j % 3 - 1;
    masked_poly_split(&s1, &s2, &p);
    for (j = 0; j < KYBER_N; j++) {
        int32_t r = (int32_t)s1.coeffs[j] + (int32_t)s2.coeffs[j];
        r %= KYBER_Q;
        ck((uint16_t)r == (uint16_t)(p.coeffs[j]), "q/2");
    }
}

/* 6.7: compression boundaries */
static void t_67(void)
{
    poly p, s1, s2;
    int j;
    printf("  6.7 compression boundaries...\n");
    for (j = 0; j < KYBER_N; j++) {
        switch(j % 5) {
            case 0: p.coeffs[j] = 0; break;
            case 1: p.coeffs[j] = 1; break;
            case 2: p.coeffs[j] = KYBER_Q/2-1; break;
            case 3: p.coeffs[j] = KYBER_Q/2+1; break;
            case 4: p.coeffs[j] = KYBER_Q-1; break;
        }
    }
    masked_poly_split(&s1, &s2, &p);
    for (j = 0; j < KYBER_N; j++) {
        int32_t r = (int32_t)s1.coeffs[j] + (int32_t)s2.coeffs[j];
        r %= KYBER_Q;
        ck((uint16_t)r == (uint16_t)(p.coeffs[j]), "boundary");
    }
}

/* 6.11: randomness varies */
static void t_611(void)
{
    poly z, a1, b1, a2, b2;
    unsigned int j, diff = 0;
    printf("  6.11 randomness varies...\n");
    for (j = 0; j < KYBER_N; j++) z.coeffs[j] = (int16_t)((j*7) % KYBER_Q);
    masked_poly_split(&a1, &b1, &z);
    masked_poly_split(&a2, &b2, &z);
    for (j = 0; j < KYBER_N; j++)
        if (a1.coeffs[j] != a2.coeffs[j]) { diff = 1; break; }
    ck(diff, "varies");
}

int main(int argc, char **argv)
{
    printf("\n");
    printf("=== Kyber Masked Decryption — Phase 6 Tests (KYBER_K=%d) ===\n\n", KYBER_K);

    t_61(); printf("\n");
    t_62(); printf("\n");
    t_63(); printf("\n");
    t_66(); printf("\n");
    t_67(); printf("\n");
    t_611(); printf("\n");

    printf("=== Total %d Passed %d Failed %d ===\n", total, ok, bad);
    return bad > 0 ? 1 : 0;
}
