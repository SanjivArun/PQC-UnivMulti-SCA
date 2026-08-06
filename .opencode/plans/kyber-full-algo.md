---
trigger: manual
---

# Kyber Full Algorithm Implementation Plan

## What is Kyber? (Explained Simply)

Kyber is a way to securely share a secret with someone over the internet. Think of it like sending a secret password through a regular mail:

- **Math magic**: Kyber uses math with big polynomials (like complex number puzzles) that are easy to do one way but nearly impossible to undo without the "key."
- **Quantum-safe**: Even quantum computers can't crack this math.
- **Standardized**: It's been officially chosen by the U.S. government as a standard (now called ML-KEM, FIPS 203).

A full Kyber operation has 3 parts:
1. **Key Generation** - Create a public key (shared openly) and a private key (kept secret)
2. **Encapsulation** - Encrypt data using the public key, producing a ciphertext and a shared secret
3. **Decapsulation** - Decrypt the ciphertext using the private key to recover the shared secret

## Current State

**Already implemented** (what KyberNTT has, ~202 lines total):
- `montgomery_reduce()` / `barrett_reduce()` — modular arithmetic
- `ntt()` / `poly_ntt()` / `polyvec_ntt()` — forward NTT
- `poly_reduce()` — Barrett reduction on polynomial coefficients

**What we need to add** — all the other math that makes Kyber work.

## Official Reference Implementation

**Repository**: [pq-crystals/kyber](https://github.com/pq-crystals/kyber) (NIST submission team)
The `ref/` directory contains the clean reference C implementation we'll copy from.

## File Changes (Minimal)

The approach is to **modify 4 existing files** and **add ~6 new files**. No files are deleted or renamed.

### Modify 4 existing files:

| File | What's Added |
|------|-------------|
| **ntt.c / ntt.h** | Add `fqmul()`, `basemul()`, `invntt()` (inverse NTT needed for decryption) |
| **poly.c / poly.h** | Add `poly_invntt_tomont()`, `poly_tomont()`, `poly_basemul_montgomery()`, `poly_add()`, `poly_sub()`, `poly_compress()`, `poly_decompress()`, `poly_tobytes()`, `poly_frombytes()`, `poly_frommsg()`, `poly_tomsg()`, `poly_getnoise_eta1()`, `poly_getnoise_eta2()`, `gen_matrix()` |
| **polyvec.c / polyvec.h** | Add `polyvec_compress()`, `polyvec_decompress()`, `polyvec_tobytes()`, `polyvec_frombytes()`, `polyvec_invntt_tomont()`, `polyvec_basemul_acc_montgomery()`, `polyvec_reduce()`, `polyvec_add()` |
| **project.c** | Rewritten as a **single sequential Kyber flow** in one `starttrigger()`/`endtrigger()` block |

### Create ~6 new files:

| Functionality | New Files | Why |
|--------------|-----------|-----|
| CBD noise sampling | **cbd.c** + **cbd.h** | Generates random noise for Kyber security |
| Constant-time helpers | **verify.c** + **verify.h** | Byte comparison and conditional copy (for CCA security) |
| SHA3/Keccak | **fips202.c** + **fips202.h** | The hash functions Kyber needs (SHA3-256, SHA3-512, SHAKE128, SHAKE256) |
| Symmetric layer | **symmetric.h** + **symmetric-shake.c** | Wraps SHA3/SHAKE for Kyber's specific needs (PRF, matrix absorb) |
| IND-CPA encryption | **indcpa.c** + **indcpa.h** | Matrix generation + actual encryption/decryption |
| KEM top-level | **kem.c** + **kem.h** | Final key generation, encapsulation, decapsulation |

### Additional new files:

| File | Purpose |
|------|---------|
| **randombytes.c** + **randombytes.h** | OS-independent CSPRNG using SHAKE-256 |

### Build file:

| File | Action |
|------|--------|
| **Makefile** | Add new .c files to SOURCES, OBJECTS (HEADERS unchanged — only `params.h` needed) |

## Implementation Order (Step-by-Step)

Functions must be built in dependency order — you can't use a function before it exists.

### Step 1: Fix ntt.c — the inverse of forward NTT

Currently has only `ntt()` (forward). Need to add:
- **`fqmul()`** — multiply two numbers mod q (wraps `montgomery_reduce`)
- **`basemul()`** — multiply two polynomials in NTT domain (uses `fqmul`)
- **`invntt()`** — inverse NTT (uses `barrett_reduce` + `fqmul` + `basemul`)

The inverse NTT is like playing a video backwards. If `ntt()` turns a polynomial into NTT form, `invntt()` converts it back.

### Step 2: Extend poly.c — polynomial operations

Currently has `poly_ntt()` and `poly_reduce()`. Need to add:
- **`poly_invntt_tomont()`** — inverse NTT on a polynomial (uses `invntt`)
- **`poly_tomont()`** — convert coefficients to Montgomery form (uses `montgomery_reduce`)
- **`poly_basemul_montgomery()`** — multiply two polynomials using NTT math (uses `basemul`)
- **`poly_add()`** / **`poly_sub()`** — add or subtract two polynomials (simple loops)
- **`poly_compress()`** / **`poly_decompress()`** — shrink/recover polynomial coefficients for sending over the wire
- **`poly_tobytes()`** / **`poly_frombytes()`** — convert a polynomial to/from bytes (for storing keys)
- **`poly_frommsg()`** / **`poly_tomsg()`** — turn a 32-byte message into a polynomial and back
- **`poly_getnoise_eta1()`** / **`poly_getnoise_eta2()`** — sample noise polynomials from seed (uses `poly_cbd` + `prf`)
- **`gen_matrix()`** — generate random matrix A deterministically from seed (uses `shake128` + `rej_uniform`)

The compression/decompression is like ZIP for math — it shrinks data so it's faster to send over a network.

### Step 3: Extend polyvec.c — vector-of-polynomials operations

Currently has only `polyvec_ntt()`. Need to add:
- **`polyvec_invntt_tomont()`** — inverse NTT on all polynomials in a vector
- **`polyvec_basemul_acc_montgomery()`** — matrix-vector multiply in NTT domain (the heart of Kyber matrix math)
- **`polyvec_reduce()`** — reduce all coefficients in all polynomials
- **`polyvec_add()`** — add two polynomial vectors
- **`polyvec_compress()`** / **`polyvec_decompress()`** — compress/decompress entire vectors
- **`polyvec_tobytes()`** / **`polyvec_frombytes()`** — serialize/deserialize vectors

A polyvec is just an array of polynomials. These functions apply the same operations to every polynomial in the array.

### Step 4: Create cbd.c — noise sampling

Centered Binomial Distribution (CBD) generates the "noise" that makes Kyber secure. Without this, anyone could break the encryption.
- **`load32_littleendian()`** / **`load24_littleendian()`** — read bytes from memory
- **`cbd2()`** / **`cbd3()`** — generate coefficients from random bits
- **`poly_cbd_eta1()`** / **`poly_cbd_eta2()`** — create polynomials from CBD

This is like rolling dice. The randomness is what keeps Kyber secure. If the dice were fixed (not random), attackers could figure out the secret key.

### Step 5: Create verify.c — constant-time helpers

These prevent "side channels" — leaking secret information by how long operations take.
- **`verify()`** — compare two byte arrays in constant time
- **`cmov()`** — copy bytes from A to B, but only if a flag is 1 (constant-time)
- **`cmov_int16()`** — same but for single integer values

Think of these like a locked door that always takes the same amount of time to open whether you have the key or not. That way someone can't tell if you tried the right key by how long it took.

### Step 6: Create fips202.c + symmetric-shake.c — hash functions

Kyber uses SHA3 and SHAKE hashes. These are cryptographic "stamp machines" that turn any input into a fixed-size unique fingerprint.

**fips202.c** — Keccak/SHA3 core:
- **`KeccakF1600_StatePermute()`** — the core permutation (24 rounds)
- **`shake128()`** / **`shake256()`** — extendable output functions (XOF)
- **`sha3_256()`** / **`sha3_512()`** — hash functions

**symmetric-shake.c** — Kyber-specific wrappers:
- **`kyber_shake128_absorb()`** — Kyber-specific SHAKE128 matrix generation
- **`kyber_shake256_prf()`** — SHAKE256 as PRF for noise generation
- **`kyber_shake256_rkprf()`** — rejection-key PRF for CCA security

### Step 7: Create indcpa.c — matrix operations + encryption

This is where all the previous pieces come together:
- **`rej_uniform()`** — rejection sampling (generate random numbers mod q)
- **`pack_pk()` / `unpack_pk()`** — serialize/deserialize public key
- **`pack_sk()` / `unpack_sk()`** — serialize/deserialize secret key
- **`pack_ciphertext()` / `unpack_ciphertext()`** — serialize/deserialize ciphertext
- **`indcpa_keypair_derand()`** — generate key pair (uses matrix A, noise, NTT multiply)
- **`indcpa_enc()`** — encrypt a message (uses A transposed, noise, NTT multiply, inverse NTT)
- **`indcpa_dec()`** — decrypt a ciphertext (uses secret key, NTT multiply, inverse NTT)

Think of the matrix A as a public puzzle everyone can see. But solving it to find the secret is nearly impossible — that's what makes Kyber secure.

### Step 8: Create kem.c — top-level KEM

This is the final "wrapper" that ties everything together:
- **`crypto_kem_keypair()`** — create keys (calls `indcpa_keypair_derand`)
- **`crypto_kem_enc()`** — encapsulate (generate random data, encrypt, produce ciphertext + shared secret)
- **`crypto_kem_dec()`** — decapsulate (decrypt ciphertext, produce shared secret, verify it's valid)

This is the public API — what other software would call to use Kyber.

### Step 9: Rewrite project.c — single sequential flow in one trigger block

All Kyber functions are called one after another inside a single:

```c
starttrigger();

// --- KEY GENERATION ---
// 1. Generate the public matrix A
// 2. Sample secret key noise
// 3. Sample error noise
// 4. Compute pk = A * secret_key + error  (NTT multiply)
// 5. Save public key and secret key

// --- ENCAPSULATION ---
// 6. Hash the public key
// 7. Derive a message and encryption randomness
// 8. Encode message as a polynomial
// 9. Sample more noise
// 10. Compute ciphertext components (NTT multiply + inverse NTT)
// 11. Compress and serialize ciphertext

// --- DECAPSULATION ---
// 12. Decompress and deserialize ciphertext
// 13. Multiply by secret key (NTT domain)
// 14. Inverse NTT to get back to normal form
// 15. Subtract from ciphertext to get message
// 16. Decode message from polynomial
// 17. Re-encrypt to verify correctness
// 18. Compare ciphertexts (constant-time)
// 19. Output final shared key

endtrigger();
```

## Sequential Function Flow in project.c

The project.c will follow the actual Kyber protocol order:

```
starttrigger();

// ============================================
// PHASE 1: KEY GENERATION (crypto_kem_keypair*)
// ============================================
// 1. randombytes → derive keygen coins
// 2. hash_g → derive public seed + noise seed
// 3. gen_matrix → generate random matrix A
// 4. poly_getnoise_eta1 × K → sample secret key polynomials
// 5. poly_getnoise_eta1 × K → sample error polynomials
// 6. polyvec_ntt × 2 → convert sk + error to NTT domain
// 7. polyvec_basemul_acc_montgomery × K → compute A * sk
// 8. poly_tomont + polyvec_add → pk = A*sk + e
// 9. polyvec_reduce → reduce pk coefficients
// 10. hash_h → hash public key for CCA security
// 11. pack_sk + pack_pk → save keys

// ============================================
// PHASE 2: ENCAPSULATION (crypto_kem_enc*)
// ============================================
// 12. randombytes → derive encryption coins
// 13. hash_h → hash public key
// 14. hash_g → derive message + encryption random coins
// 15. poly_frommsg → encode message as polynomial
// 16. poly_getnoise × (3K+1) → sample encryption randomness
// 17. polyvec_ntt → convert s to NTT domain
// 18. polyvec_basemul_acc_montgomery → compute b = A^T * s + ep
// 19. polyvec_basemul_acc_montgomery → compute v = pk * s + epp + msg
// 20. polyvec_invntt_tomont + poly_invntt_tomont → convert back
// 21. poly_add + polyvec_reduce → clean up
// 22. polyvec_compress + poly_compress → compress ciphertext
// 23. pack_ciphertext → serialize ct

// ============================================
// PHASE 3: DECAPSULATION (crypto_kem_dec)
// ============================================
// 24. unpack_ciphertext → decompress b, v
// 25. polyvec_ntt → convert b to NTT domain
// 26. polyvec_basemul_acc_montgomery → compute mp = sk * b
// 27. poly_invntt_tomont → convert back
// 28. poly_sub → diff = v - mp
// 29. poly_reduce → reduce coefficients
// 30. poly_tomsg → decode to 32-byte message
// 31. hash_g → re-derive encryption coins
// 32. indcpa_enc → re-encrypt to verify
// 33. verify → constant-time comparison
// 34. rkprf → generate rejection key
// 35. cmov → select correct shared secret
```

## What Each Phase Does (Simple Explanation)

```
Phase 1: KEY GENERATION
  A          sk              e
  ┌───┐     ┌───┐           ┌───┐
  │ 1 │     │ 0 │           │ 3 │   → random "noise" that makes crypto secure
  └───┘     └───┘           └───┘
  ┌───┐     ┌───┐           ┌───┐
  │ 2 │     │ 0 │           │ 3 │
  └───┘     └───┘           └───┘

  Public Key = A × sk + e    →  this is safe to share with everyone
  Secret Key = sk             →  this stays private

Phase 2: ENCAPSULATION (encrypting)
  Public Key + Message + Randomness
  ═══════════════════════════════════
  Ciphertext = A^T × s + ep       (noise makes it look random)
  SharedSecret = v + epp + msg    (this becomes the shared encryption key)

Phase 3: DECAPSULATION (decrypting)
  Secret Key + Ciphertext
  ═══════════════════════
  Recovered Message = v - (sk × C₁) + msg
                    = msg                 (noise cancels out!)
  Check: Re-encrypt and compare to verify correctness.
```

## Files Summary

**Modified (4 files):**
- `ntt.c` (+ `ntt.h`) — add fqmul, basemul, invntt
- `poly.c` (+ `poly.h`) — add ~15 new poly functions
- `polyvec.c` (+ `polyvec.h`) — add 8 new polyvec functions
- `project.c` — complete rewrite with full flow in single trigger

**Created (~8 files):**
- `cbd.c` + `cbd.h` — noise sampling (~30 lines)
- `verify.c` + `verify.h` — constant-time utilities (~25 lines)
- `fips202.c` + `fips202.h` — SHA3/Keccak (~400 lines from reference)
- `symmetric.h` + `symmetric-shake.c` — symmetric wrappers (~50 lines)
- `indcpa.c` + `indcpa.h` — IND-CPA encryption (~200 lines)
- `kem.c` + `kem.h` — KEM wrapper (~130 lines)
- `randombytes.c` + `randombytes.h` — deterministic CSPRNG (~20 lines)

**Total:** ~14 new files + 4 modified files + Makefile update

## Important Numbers

- **KYBER_K = 2** (Kyber-512 security level — smallest/fastest)
- **q = 3329** (the modulus — a prime number)
- **n = 256** (each polynomial has 256 coefficients)
- **eta1 = 3** (noise parameter for Kyber-512)
- **eta2 = 2** (always 2)
- **Shared secret = 32 bytes** (256 bits)
- **Public key ≈ 800 bytes**
- **Secret key ≈ 1,632 bytes**
- **Ciphertext ≈ 768 bytes**
