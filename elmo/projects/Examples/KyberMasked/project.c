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
  uint8_t kgbuf[2 * KYBER_SYMBYTES];
  uint8_t enc_coins[KYBER_SYMBYTES];
  uint8_t enc_kr[2 * KYBER_SYMBYTES];

  read2bytes(&nb_challenges);

  for(num_challenge = 0; num_challenge < nb_challenges; num_challenge++) {

    // Prepare keygen input
    for(i = 0; i < KYBER_SYMBYTES; i++)
      kgbuf[i] = num_challenge + i;
    for(i = 0; i < KYBER_SYMBYTES; i++)
      kgbuf[KYBER_SYMBYTES + i] = 0xAA;

    // Prepare encap input
    for(i = 0; i < KYBER_SYMBYTES; i++)
      enc_coins[i] = i + 0x55;

    // --- KEY GENERATION ---
    starttrigger();
    hash_g(kgbuf, kgbuf, 2 * KYBER_SYMBYTES);
    crypto_kem_keypair_derand(pk, sk, kgbuf);
    endtrigger();

    // Print public key and secret key
    for(i = 0; i < KYBER_PUBLICKEYBYTES; i++) {
      u16 = (uint16_t)pk[i];
      print2bytes(&u16);
    }
    for(i = 0; i < KYBER_SECRETKEYBYTES; i++) {
      u16 = (uint16_t)sk[i];
      print2bytes(&u16);
    }

    // --- ENCAPSULATION ---
    starttrigger();
    crypto_kem_enc_derand(ct, enc_kr, pk, enc_coins);
    endtrigger();

    // Print shared secret and ciphertext
    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_a[i];
      print2bytes(&u16);
    }
    for(i = 0; i < KYBER_CIPHERTEXTBYTES; i++) {
      u16 = (uint16_t)ct[i];
      print2bytes(&u16);
    }

    // --- DECAPSULATION ---
    starttrigger();
    crypto_kem_dec(ss_b, ct, sk);
    endtrigger();

    // Print decrypted shared secret
    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_b[i];
      print2bytes(&u16);
    }
  }

  endprogram();

  return 0;
}
