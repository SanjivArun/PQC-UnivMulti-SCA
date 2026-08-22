#include <stdio.h>
#include <stdlib.h>

#include "elmoasmfunctionsdef-extension.h"
#include "polyvec.h"
#include "params.h"

int main(void) {
  uint16_t num_challenge, nb_challenges;
  int j, k;
  polyvec skpv;

  read2bytes(&nb_challenges);
  for(num_challenge = 0; num_challenge < nb_challenges; num_challenge++) {
    // consume the entire challenge exactly as projectclass.py writes it
    // (KYBER_K * KYBER_N 16-bit values), same pattern as the working KyberNTT
    for(j = 0; j < KYBER_K; j++)
      for(k = 0; k < KYBER_N; k++)
        read2bytes((uint16_t*)&skpv.vec[j].coeffs[k]);

    // leaking operation wrapped in the trigger window
    starttrigger();
    polyvec_ntt(&skpv);
    endtrigger();

    for(j = 0; j < KYBER_K; j++)
      for(k = 0; k < KYBER_N; k++)
        print2bytes((uint16_t*)&skpv.vec[j].coeffs[k]);
  }

  endprogram();
  return 0;
}
