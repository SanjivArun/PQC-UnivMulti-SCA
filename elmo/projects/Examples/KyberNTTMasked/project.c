#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "elmoasmfunctionsdef-extension.h"

// ELMO API :
//  - printbyte(addr): Print single byte located at address 'addr' to output file;
//  - randbyte(addr): Load byte of random to memory address 'addr';
//  - readbyte(addr): Read byte from input file to address 'addr'.
// ELMO API (extension) :
//  - print2bytes, rand2bytes and read2bytes: idem, but for an address pointing on 2 bytes;
//  - print4bytes, rand4bytes and read4bytes: idem, but for an address pointing on 4 bytes.

#include "polyvec.h"
#include "params.h"
#include "reduce.h"

int main(void) {
  uint16_t num_challenge, nb_challenges;
  int j, k;
  polyvec skpv, s0, s1;

  read2bytes(&nb_challenges);
  for(num_challenge=0; num_challenge<nb_challenges; num_challenge++) {

    // Load the private vector s
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++)
        read2bytes((uint16_t*) &skpv.vec[j].coeffs[k]);

    // MASKED VERSION (the original unmasked body is kept below in comments):
    //
    // Original code:
    //   starttrigger(); // To start a new trace
    //   // Do the leaking operations here...
    //   polyvec_ntt(&skpv);
    //   endtrigger(); // To end the current trace
    //
    // Masking: split s into two shares s0, s1 with s = s0 + s1 (mod q).
    // s0 is random, s1 is derived. The split happens BEFORE starttrigger()
    // so that the true secret s never appears inside the recorded trace.
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++) {
        rand2bytes((uint16_t*) &s0.vec[j].coeffs[k]);
        s0.vec[j].coeffs[k] = barrett_reduce(s0.vec[j].coeffs[k]); /* s0 in [0, q] */
        s1.vec[j].coeffs[k] = (int16_t)(skpv.vec[j].coeffs[k] - s0.vec[j].coeffs[k]);
        while(s1.vec[j].coeffs[k] < 0)
          s1.vec[j].coeffs[k] += KYBER_Q; /* s1 in [0, q] */
      }

    starttrigger(); // To start a new trace

    // Leaking operations... now only on the random-looking shares.
    // NTT is linear over Z_q: NTT(s) = NTT(s0) + NTT(s1).
    polyvec_ntt(&s0);
    polyvec_ntt(&s1);

    endtrigger(); // To end the current trace

    // Recombine the two shares (mod q) OUTSIDE the recorded trace:
    // writing the true NTT(s) inside starttrigger()/endtrigger() would leak it.
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++)
        skpv.vec[j].coeffs[k] = barrett_reduce(
            (int16_t)(s0.vec[j].coeffs[k] + s1.vec[j].coeffs[k]));

    // Print the results of the computation
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++)
        print2bytes((uint16_t*) &skpv.vec[j].coeffs[k]);
  }

  endprogram(); // To indicate to ELMO that the simulation is finished

  return 0;
}
