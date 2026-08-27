/*
 * Kyber Masked Decryption - Phase 4: Masked Compression (A2B)
 *
 * masked_poly_tomsg() computes Compress_q(x, 1) from two ARITHMETIC shares
 * x1, x2 (x1 + x2 = x mod q) WITHOUT ever reconstructing x = x1 + x2.
 *
 * Reference: the masked compaction from "First-Order Masked Kyber on ARM
 * Cortex-M4" (HKK22, ePrint 2022/058), which is the mkm4 implementation.
 * It uses:
 *   - the OSFP18 mod-q -> mod-2^16 transform (no rejection sampling), and
 *   - the BDV21 "Debrāize (fixed)" single-lookup table-based
 *     arithmetic-to-Boolean conversion (ePrint 2021/067, Algorithm 8),
 *   which is the secure variant (independent per-iteration mask r_i).
 *
 * Configuration (mkm4): k = 8, a2bn = 2 (two 8-bit iterations, 16-bit
 * value), mod2k = 0xff, mod2nk = 0xffff. The table has a2bn * 2^(k+1)
 * = 1024 16-bit entries.
 */

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "masked_a2b.h"
#include "randombytes.h"

#define A2B_K       8        /* bits per chunk          */
#define A2B_N       2        /* number of chunks        */
#define A2B_MOD2K   0xffu    /* mask for one k-bit chunk */
#define A2B_MOD2NK  0xffffu  /* mask for the full value  */

/* Flat index into the precomputed table (mkm4 def.h index()). */
#define A2B_INDEX(i, beta, A) \
    (((uint16_t)(i) << 9) | ((uint16_t)(beta) << A2B_K) | (uint16_t)(A))

/* ------------------------------------------------------------------ */
/*  Fresh 32-bit random (advances the deterministic PRNG, constant
 *  cycle count, no rejection sampling -> constant ELMO trace length). */
/* ------------------------------------------------------------------ */
static uint32_t randomint(void)
{
    uint8_t b[4];
    uint32_t v;

    randombytes(b, 4);
    v = ((uint32_t)b[0] << 24)
      | ((uint32_t)b[1] << 16)
      | ((uint16_t)b[2] << 8)
      | (uint32_t)b[3];
    return v;
}

/* ------------------------------------------------------------------ */
/*  BDV21 single-lookup A2B context (mkm4 a2b_singlelookup.c)          */
/* ------------------------------------------------------------------ */
static uint16_t A2B_T[A2B_N * 512];          /* 1024 entries = 2 KiB   */

typedef struct {
    uint16_t *T;
    uint8_t  rho;
    uint16_t rrr;
} a2b_ctx_t;

static a2b_ctx_t A2B_ctx;

static void A2B_init(void);
static uint32_t A2B_convert_sw(uint32_t A, uint32_t R);

/* Table generation (BDV21 Algorithm 7 — "Debrāize (fixed)").
 * Independent per-iteration mask r_i; single shared carry mask rho. */
static void A2B_init(void)
{
    uint8_t r[A2B_N];
    uint32_t buff;
    uint32_t rrr = 0;
    size_t i;
    size_t A;

    A2B_ctx.rho = (uint8_t)(randomint() & 1u);
    buff = randomint() & 0xFFFFu;

    for (i = 0; i < A2B_N; i++) {
        r[i] = (uint8_t)((buff >> (i * A2B_K)) & A2B_MOD2K);
        rrr |= ((uint16_t)r[i] << (i * A2B_K));
    }
    A2B_ctx.rrr = (uint16_t)rrr;

    for (i = 0; i < A2B_N; i++) {
        for (A = 0; A < 256; A++) {
            A2B_T[A2B_INDEX(i, A2B_ctx.rho, A)] =
                (uint16_t)(A + r[i]) ^ (uint16_t)((A2B_ctx.rho << A2B_K) | r[i]);
            A2B_T[A2B_INDEX(i, (uint8_t)(A2B_ctx.rho ^ 1), A)] =
                (uint16_t)(A + r[i] + 1) ^ (uint16_t)((A2B_ctx.rho << A2B_K) | r[i]);
        }
    }
}

/* Convert arithmetic (A, R) with A + R = x mod 2^16 into a Boolean
 * mask B with B XOR rrr = x. Never sums A and R (Algorithm 8). */
static uint32_t A2B_convert_sw(uint32_t A, uint32_t R)
{
    size_t i;
    uint32_t A_l;
    uint32_t R_l;
    uint32_t betaBi;
    uint32_t beta = A2B_ctx.rho;
    uint32_t Bi;
    uint32_t B = 0;

    A = (A - A2B_ctx.rrr) & A2B_MOD2NK;

    for (i = 0; i < A2B_N; i++) {
        uint32_t mask = (uint32_t)((1u << (A2B_N - i) * A2B_K) - 1);
        R_l = R & A2B_MOD2K;
        A = (A + R_l) & mask;
        A_l = A & A2B_MOD2K;

        betaBi = A2B_T[A2B_INDEX(i, (uint32_t)beta, A_l)];
        Bi     = betaBi & A2B_MOD2K;
        beta   = betaBi >> A2B_K;

        Bi = Bi ^ R_l;
        B  |= (Bi << (i * A2B_K));

        A >>= A2B_K;
        R >>= A2B_K;
    }

    return B ^ A2B_ctx.rrr;
}

static uint32_t A2B_convert(uint32_t A, uint32_t R)
{
    return A2B_convert_sw(A, R);
}

/* ------------------------------------------------------------------ */
/*  Masked poly_tomsg (mkm4 crypto_kem/kyber768/m4/masked-poly.c)      */
/* ------------------------------------------------------------------ */
/*
 * For each coefficient of the two arithmetic shares (a->coeffs[0]=a1,
 * a->coeffs[1]=a2), a1 + a2 = x mod q. The bit Compress_q(x, 1) is
 * recovered from the Boolean shares WITHOUT reconstructing x:
 *
 *   1. a1 <- a1 - q/4 (mod q)     (threshold shift on share 0 only)
 *   2. OSFP18 transform: mod-q shares -> mod-2^16 shares (no rejection)
 *   3. c0  <- c0 - q/2           (threshold shift on share 0 only)
 *   4. A2B of the mod-2^16 share, take the MSB (bit 15) of each share
 *
 * Message bit = (share0 bit) XOR (share1 bit) = Compress_q(x, 1).
 */
void masked_poly_tomsg(uint8_t m1[KYBER_INDCPA_MSGBYTES],
                       uint8_t m2[KYBER_INDCPA_MSGBYTES],
                       const poly *x1,
                       const poly *x2)
{
    size_t i;         /* symbol index (0 .. MSGBYTES-1) */
    size_t j;         /* bit index within a byte (0 .. 7) */
    uint16_t a1;
    uint16_t a2;
    uint16_t random;  /* unsigned 16-bit mask (0 .. 2^16-1) */
    uint16_t k11;
    uint16_t k12;
    uint16_t k21;
    uint16_t k22;
    uint16_t y1;
    uint16_t y2;
    uint16_t z0;
    uint16_t z[2];
    uint16_t c0;
    uint16_t c1;
    uint32_t tmp;

    for (i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
        m1[i] = 0;
        m2[i] = 0;
    }

    for (i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
        for (j = 0; j < 8; j++) {

            A2B_init();          /* fresh table + rrr per coefficient */

            /* Constant-time mod-q reduction (no division): the value
             * x - q/4 + q lies in [-q/4, 3q/2], so one conditional
             * subtraction yields the result in [0, q). Using % would
             * call __aeabi_idivmod, whose iteration count depends on
             * the dividend, making the ELMO trace length data-dependent. */
            int32_t a1v = (int32_t)x1->coeffs[8 * i + j] - KYBER_Q / 4 + KYBER_Q;
            if (a1v >= KYBER_Q) a1v -= KYBER_Q;
            if (a1v < 0) a1v += KYBER_Q;
            a1 = (uint16_t)a1v;
            a2 = (uint16_t)((int32_t)x2->coeffs[8 * i + j]);

            /* --- OSFP18 mod-q -> mod-2^16 transform --- */
            tmp    = randomint();
            y1     = (uint16_t)(tmp & 0xFFFFu);
            random = (uint16_t)((tmp >> 16) & 0xFFFFu);

            y2 = (uint16_t)(((uint32_t)a1 + (0xFFFFu - y1) + 1));
            y2 = (uint16_t)(y2 + a2);
            z0 = (uint16_t)(((uint32_t)y1 + (0xFFFFu - KYBER_Q) + 1));

            z[0] = A2B_convert(z0, y2);
            z[1] = y2;

            tmp = randomint();
            k11 = (uint16_t)(tmp & 0xFFFFu);
            k21 = (uint16_t)((tmp >> 16) & 0xFFFFu);

            k12 = (uint16_t)(((uint32_t)((uint16_t)z[0] >> 15 ^ 1u) + (0xFFFFu - k11) + 1));
            k22 = (uint16_t)(((uint32_t)((uint16_t)z[1] >> 15)     + (0xFFFFu - k21) + 1));

            c0 = (uint16_t)(((uint32_t)random + y1));
            c0 = (uint16_t)((uint32_t)c0
                   - (uint32_t)(((uint16_t)z[0] >> 15) ^ 1) * (uint32_t)KYBER_Q);
            c1 = (uint16_t)((uint32_t)y2 + (0xFFFFu - random) + 1);
            c0 = (uint16_t)((uint32_t)c0
                   - (uint32_t)((uint16_t)z[1] >> 15) * (uint32_t)KYBER_Q);
            /* Each 2*k*u*v*Q term kept mod 2^16 with 32-bit intermediates so no
             * slow 64-bit multiplier is ever needed on the Cortex-M0 target. */
            { uint16_t t = (uint16_t)(((uint32_t)((uint16_t)2 * KYBER_Q) * k11));
              t = (uint16_t)(((uint32_t)t * k21));
              c0 = (uint16_t)((uint32_t)c0 + t); }
            { uint16_t t = (uint16_t)(((uint32_t)((uint16_t)2 * KYBER_Q) * k11));
              t = (uint16_t)(((uint32_t)t * k22));
              c0 = (uint16_t)((uint32_t)c0 + t); }
            { uint16_t t = (uint16_t)(((uint32_t)((uint16_t)2 * KYBER_Q) * k12));
              t = (uint16_t)(((uint32_t)t * k21));
              c0 = (uint16_t)((uint32_t)c0 + t); }
            { uint16_t t = (uint16_t)(((uint32_t)((uint16_t)2 * KYBER_Q) * k12));
              t = (uint16_t)(((uint32_t)t * k22));
              c0 = (uint16_t)((uint32_t)c0 + t); }
            /* --- end OSFP18 transform --- */

            c0 = (uint16_t)((uint32_t)c0 + (0xFFFFu - KYBER_Q / 2) + 1);

            c0 = A2B_convert(c0, c1);

            m1[i] += (uint8_t)(((c0 >> 15) & 1) << j);
            m2[i] += (uint8_t)(((c1 >> 15) & 1) << j);
        }
    }
}
