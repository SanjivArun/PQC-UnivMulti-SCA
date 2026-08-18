/* Deterministic PRNG for testing — replaces randombytes.c in build */
#include <stdint.h>
#include "randombytes.h"

static uint64_t prng_a = 6364136223846793005ULL;
static uint64_t prng_c = 1442695040888963407ULL;
static uint64_t prng_x = 42;

void randombytes(uint8_t *out, size_t outlen)
{
    size_t i;
    for (i = 0; i < outlen; i++) {
        prng_x = prng_x * prng_a + prng_c;
        out[i] = (uint8_t)(prng_x >> 48);
    }
}

void randombytes_seed(uint64_t seed)
{
    prng_x = seed;
}
