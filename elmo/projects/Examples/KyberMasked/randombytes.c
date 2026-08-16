#include <stdint.h>
#include <stdlib.h>
#include "fips202.h"
#include "randombytes.h"

void randombytes(uint8_t *out, size_t outlen)
{
  keccak_state s;
  size_t nblocks;

  shake256_absorb_once(&s, (uint8_t*)&outlen, sizeof(outlen));

  nblocks = outlen / SHAKE256_RATE;
  shake256_squeezeblocks(out, nblocks, &s);
  out += nblocks * SHAKE256_RATE;
  outlen -= nblocks * SHAKE256_RATE;

  shake256_squeeze(out, outlen, &s);
}
