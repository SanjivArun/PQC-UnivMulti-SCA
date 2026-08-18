#include <stdio.h>
#include <stdlib.h>

#include "elmoasmfunctionsdef-extension.h"
#include "kem.h"
#include "symmetric.h"

int main(void) {
  uint8_t pk[KYBER_PUBLICKEYBYTES];
  uint8_t sk[KYBER_SECRETKEYBYTES];
  uint8_t ct[KYBER_CIPHERTEXTBYTES];
  uint8_t ss_a[KYBER_SSBYTES];
  uint8_t ss_b[KYBER_SSBYTES];
  uint16_t num_challenge, nb_challenges;
  uint8_t i;
  uint16_t u16;

  read2bytes(&nb_challenges);

  for(num_challenge = 0; num_challenge < nb_challenges; num_challenge++) {

    crypto_kem_enc(ct, ss_a, pk);
    crypto_kem_dec(ss_b, ct, sk);

    starttrigger();

    // --- KEY GENERATION ---
    uint8_t kgbuf[2 * KYBER_SYMBYTES];
    for(i = 0; i < KYBER_SYMBYTES; i++)
      kgbuf[i] = num_challenge + i;
    for(i = 0; i < KYBER_SYMBYTES; i++)
      kgbuf[KYBER_SYMBYTES + i] = 0xAA;

    hash_g(kgbuf, kgbuf, 2 * KYBER_SYMBYTES);
    crypto_kem_keypair_derand(pk, sk, kgbuf);

    // --- ENCAPSULATION ---
    uint8_t enc_coins[KYBER_SYMBYTES];
    uint8_t enc_kr[2 * KYBER_SYMBYTES];
    for(i = 0; i < KYBER_SYMBYTES; i++)
      enc_coins[i] = i + 0x55;

    crypto_kem_enc_derand(ct, enc_kr, pk, enc_coins);

    // --- DECAPSULATION ---
    crypto_kem_dec(ss_b, ct, sk);

    endtrigger();

    // Print shared secrets
    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_a[i];
      print2bytes(&u16);
    }

    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_b[i];
      print2bytes(&u16);
    }

    // Print ciphertext (16-bit values)
    u16 = 0;
    for(num_challenge = 0; num_challenge < KYBER_CIPHERTEXTBYTES; num_challenge++) {
      u16 = (uint16_t)ct[num_challenge];
      print2bytes(&u16);
    }
  }

  endprogram();

  return 0;
}
