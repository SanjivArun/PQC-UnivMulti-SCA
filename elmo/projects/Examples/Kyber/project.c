#include <stdio.h>
#include <stdlib.h>

#include "elmoasmfunctionsdef-extension.h"
#include "kem.h"
#include "indcpa.h"
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

    // ============================================================
    // --- KEY GENERATION ---
    // crypto_kem_keypair_derand(pk, sk, kgbuf)
    // ============================================================

    // hash_g(kgbuf, kgbuf, 2 * KYBER_SYMBYTES)
    // starttrigger();
    hash_g(kgbuf, kgbuf, 2 * KYBER_SYMBYTES);
    // endtrigger();

    // indcpa_keypair_derand(pk, sk, kgbuf)
    // starttrigger();
    indcpa_keypair_derand(pk, sk, kgbuf);
    // endtrigger();

    // memcpy(sk+KYBER_INDCPA_SECRETKEYBYTES, pk, KYBER_PUBLICKEYBYTES)
    // starttrigger();
    for(i = 0; i < KYBER_PUBLICKEYBYTES; i++)
      sk[KYBER_INDCPA_SECRETKEYBYTES + i] = pk[i];
    // endtrigger();

    // hash_h(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES)
    // starttrigger();
    hash_h(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES);
    // endtrigger();

    // memcpy(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, coins+KYBER_SYMBYTES, KYBER_SYMBYTES)
    // starttrigger();
    for(i = 0; i < KYBER_SYMBYTES; i++)
      sk[KYBER_SECRETKEYBYTES-KYBER_SYMBYTES + i] = kgbuf[KYBER_SYMBYTES + i];
    // endtrigger();

    // Print public key and secret key
    for(i = 0; i < KYBER_PUBLICKEYBYTES; i++) {
      u16 = (uint16_t)pk[i];
      print2bytes(&u16);
    }
    for(i = 0; i < KYBER_SECRETKEYBYTES; i++) {
      u16 = (uint16_t)sk[i];
      print2bytes(&u16);
    }

    // ============================================================
    // --- ENCAPSULATION ---
    // crypto_kem_enc_derand(ct, enc_kr, pk, enc_coins)
    // ============================================================

    // hash_h(buf+KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES)
    // hash_g(kr, buf, 2*KYBER_SYMBYTES)
    // indcpa_enc(ct, buf, pk, kr+KYBER_SYMBYTES)
    // memcpy(ss,kr,KYBER_SYMBYTES)
    // starttrigger();
    crypto_kem_enc_derand(ct, enc_kr, pk, enc_coins);
    // endtrigger();

    // Print shared secret and ciphertext
    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_a[i];
      print2bytes(&u16);
    }
    for(i = 0; i < KYBER_CIPHERTEXTBYTES; i++) {
      u16 = (uint16_t)ct[i];
      print2bytes(&u16);
    }

    // ============================================================
    // --- DECAPSULATION ---
    // crypto_kem_dec(ss_b, ct, sk)
    // ============================================================

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
