#include <stdio.h>
#include <stdlib.h>

#include "elmoasmfunctionsdef-extension.h"
#include "masked_indcpa.h"
#include "params.h"
#include "polyvec.h"
#include "poly.h"

/* Decryption-only TVLA (first-order masked).
 *
 * Same interface as the unmasked project, but the trigger window wraps the
 * first-order masked INDCPA decapsulation (masked_indcpa_dec): the secret s is
 * split into arithmetic shares, the dot product is computed on the shares, and
 * the message decode uses A2B. A sound first-order mask should make this
 * version show ~0 leakage, whereas the unmasked version shows a large count. */

int main(void) {
  uint8_t skbytes[KYBER_INDCPA_SECRETKEYBYTES];
  uint8_t ctbytes[KYBER_INDCPA_BYTES];
  uint8_t mbuf[KYBER_INDCPA_MSGBYTES];
  uint16_t nb_challenges;
  uint16_t num_challenge;
  uint16_t j;
  uint16_t k;
  uint16_t i;
  uint16_t u16;
  static polyvec s;

  read2bytes(&nb_challenges);

  for(num_challenge = 0; num_challenge < nb_challenges; num_challenge++) {

    /* 1. Read the challenge exactly as projectclass.py writes it. */
    for(j = 0; j < KYBER_K; j++)
      for(k = 0; k < KYBER_N; k++)
        read2bytes((uint16_t*)&s.vec[j].coeffs[k]);

    /* 2. Deterministic, non-degenerate ciphertext derived from the secret. */
    for(i = 0; i < KYBER_INDCPA_BYTES; i++)
      ctbytes[i] = (uint8_t)(0x20 + i + s.vec[i % KYBER_K].coeffs[i % KYBER_N]);

    /* 3. Secret key = NTT-domain s packed as keygen does. */
    polyvec_ntt(&s);
    polyvec_tobytes(skbytes, &s);

    /* 4. Masked decryption -- first-order leakage should be hidden here. */
    starttrigger();
    masked_indcpa_dec(mbuf, ctbytes, skbytes);
    endtrigger();

    /* 5. Side-effect so the compiler cannot drop the decryption. */
    for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
      u16 = (uint16_t)mbuf[i];
      print2bytes(&u16);
    }
  }

  endprogram();
  return 0;
}
