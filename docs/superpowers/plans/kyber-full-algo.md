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

**Already implemented** (what KyberNTT has):
- `montgomery_reduce()` - Montgomery modular reduction
- `barrett_reduce()` - Barrett modular reduction  
- `ntt()` - Forward Number Theoretic Transform
- `poly_ntt()` - Apply NTT to a polynomial
- `poly_reduce()` - Reduce polynomial coefficients
- `polyvec_ntt()` - Apply NTT to a vector of polynomials

**Missing** (what we need to add):
- Inverse NTT and all polynomial multiplication
- Compression/decompression (for sending data over the wire)
- Noise sampling (generating random "noise" for security)
- Hash functions (SHA3/SHAKE)
- Key generation, encryption, and decryption
- The full KEM (encapsulation/decapsulation)

## Official Reference Implementation

**Repository**: [pq-crystals/kyber](https://github.com/pq-crystals/kyber) (NIST submission)
The `ref/` directory contains the clean reference C implementation we'll follow.

Also see [pq-code-package/mlkem-native](https://github.com/pq-code-package/mlkem-native) for the newer standardized version.

## Implementation Order (Step-by-Step)

We need to build the functions in dependency order - you can't use a function before you create it. Here is the exact order:

### Phase 1: Polynomial Multiplication (NTT Domain)

1. **`fqmul()`** ✅ - Field multiplication using Montgomery reduction (ntt.c)
2. **`basemul()`** ✅ - Base multiplication in NTT domain (ntt.c)
3. **`invntt()`** ✅ - Inverse NTT (ntt.c)
4. **`poly_invntt_tomont()`** ✅ - Apply inverse NTT + convert to Montgomery domain (poly.c)
5. **`poly_tomont()`** ✅ - Convert coefficients to Montgomery form (poly.c)
6. **`poly_basemul_montgomery()`** ✅ - Multiply two polynomials in NTT domain (poly.c)
7. **`polyvec_basemul_acc_montgomery()`** ✅ - Dot product of two polyvecs (polyvec.c)
8. **`polyvec_invntt_tomont()`** ✅ - Apply inverse NTT to all polyvecs (polyvec.c)

### Phase 2: Polynomial Arithmetic

9. **`poly_add()`** ✅ - Add two polynomials (poly.c)
10. **`poly_sub()`** ✅ - Subtract two polynomials (poly.c)
11. **`poly_compress()`** ✅ - Compress polynomial for serialization (poly.c)
12. **`poly_decompress()`** ✅ - Decompress polynomial back (poly.c)
13. **`poly_tobytes()`** ✅ - Serialize polynomial to bytes (poly.c)
14. **`poly_frombytes()`** ✅ - Deserialize polynomial from bytes (poly.c)
15. **`poly_frommsg()`** ✅ - Encode a 32-byte message as a polynomial (poly.c)
16. **`poly_tomsg()`** ✅ - Decode polynomial back to 32-byte message (poly.c)

### Phase 3: Polyvec Compression

17. **`polyvec_compress()`** ✅ - Vector compression (polyvec.c)
18. **`polyvec_decompress()`** ✅ - Vector decompression (polyvec.c)
19. **`polyvec_tobytes()`** ✅ - Vector serialization (polyvec.c)
20. **`polyvec_frombytes()`** ✅ - Vector deserialization (polyvec.c)
21. **`polyvec_reduce()`** ✅ - Reduce all coefficients in a vector (polyvec.c)
22. **`polyvec_add()`** ✅ - Add two polyvecs (polyvec.c)

### Phase 4: Noise Sampling (CBD) + Verify

23. **`load32_littleendian()`** ✅ - Load 4 bytes as 32-bit integer (cbd.c)
24. **`load24_littleendian()`** ✅ - Load 3 bytes as 32-bit integer (cbd.c)
25. **`cbd2()`** ✅ - Centered Binomial Distribution with eta=2 (cbd.c)
26. **`cbd3()`** ✅ - CBD with eta=3 (cbd.c)
27. **`poly_cbd_eta1()`** ✅ - CBD sampling with eta1 (cbd.c)
28. **`poly_cbd_eta2()`** ✅ - CBD sampling with eta2 (cbd.c)

### Phase 5: Hash Functions (Keccak/SHA3) ✅

30. **`KeccakF1600_StatePermute()`** ✅ - Core Keccak permutation (24 rounds) (fips202.c)
31. **`shake128_absorb_once()` / `shake128_squeezeblocks()`** ✅ - SHAKE128 (fips202.c)
32. **`shake256_absorb_once()` / `shake256_squeezeblocks()`** ✅ - SHAKE256 (fips202.c)
33. **`sha3_256()`** ✅ - SHA3-256 hash (fips202.c)
34. **`sha3_512()`** ✅ - SHA3-512 hash (fips202.c)
35. **`kyber_shake128_absorb()`** ✅ - Kyber-specific SHAKE128 (symmetric-shake.c)
36. **`kyber_shake256_prf()`** ✅ - SHAKE256 as PRF for noise (symmetric-shake.c)

### Phase 6: Matrix Generation (IND-CPA) ✅

37. **`rej_uniform()`** ✅ - Rejection sampling for uniform values (indcpa.c)
38. **`gen_matrix()`** ✅ - Generate the public matrix A from a seed (indcpa.c)
39. **`pack_pk()` / `unpack_pk()`** ✅ - Serialize/deserialize public key (indcpa.c)
40. **`pack_sk()` / `unpack_sk()`** ✅ - Serialize/deserialize secret key (indcpa.c)
41. **`pack_ciphertext()` / `unpack_ciphertext()`** ✅ - Serialize/deserialize ciphertext (indcpa.c)

### Phase 7: IND-CPA Encryption ✅

42. **`indcpa_keypair_derand()`** ✅ - Generate public/secret key pair (indcpa.c)
43. **`indcpa_enc()`** ✅ - Encrypt a message (indcpa.c)
44. **`indcpa_dec()`** ✅ - Decrypt a ciphertext (indcpa.c)

### Phase 8: Constant-Time Utilities ✅

45. **`verify()`** ✅ - Constant-time byte comparison (verify.c)
46. **`cmov()`** ✅ - Constant-time conditional copy (verify.c)

### Phase 9: KEM Layer (Top-Level) ✅

47. **`crypto_kem_keypair_derand()`** ✅ - Key generation wrapper (kem.c)
48. **`crypto_kem_enc_derand()`** ✅ - Encapsulation wrapper (kem.c)
49. **`crypto_kem_dec()`** ✅ - Decapsulation wrapper (kem.c)

### Phase 10: Update project.c ✅

50. **Update `project.c`** ✅ - Call all functions in sequential order (build succeeds clean)

---

## Why These Phases? (Think of It Like Cooking)

```
Phase 1   → Chop your ingredients (basic math tools)
Phase 2   → Prepare the sauces (polynomial operations)
Phase 3   → Prepare the vegetables (vector operations)
Phase 4   → Get the spices ready (random noise)
Phase 5   → Preheat the oven (hash functions)
Phase 6   → Mix the batter (matrix operations)
Phase 7   → Bake the cake (encryption)
Phase 8   → Check it's not burnt (constant-time checks)
Phase 9   → Serve the cake (KEM wrapper)
Phase 10  → Eat it (sequential test in project.c)
```

## Sequential Function Flow in project.c

The project.c will follow the actual Kyber protocol order:

```
1. KEY GENERATION
   ├─ randombytes (get random numbers)
   ├─ hash_G (SHA3-512) → derive public seed + noise seed
   ├─ gen_matrix (create random matrix A)
   ├─ poly_getnoise × K (sample secret key noise)
   ├─ poly_getnoise × K (sample error noise)
   ├─ polyvec_ntt (convert to NTT domain)
   ├─ polyvec_basemul_acc_montgomery (multiply matrices)
   └─ polyvec_reduce (clean up results)
   └─ pack_pk, pack_sk (save keys)

2. ENCAPSULATION
   ├─ hash_H (SHA3-256 on public key)
   ├─ hash_G (derive message + encryption coins)
   ├─ poly_frommsg (encode message)
   ├─ poly_getnoise × (3K+1) (sample randomness)
   ├─ polyvec_ntt (convert to NTT domain)
   ├─ polyvec_basemul_acc_montgomery (matrix-vector multiply)
   ├─ polyvec_invntt_tomont (convert back to normal)
   ├─ poly_invntt_tomont (convert back to normal)
   ├─ poly_add, polyvec_reduce (clean up)
   └─ pack_ciphertext (compress + serialize)

3. DECAPSULATION
   ├─ unpack_ciphertext (decompress ciphertext)
   ├─ polyvec_ntt (convert to NTT domain)
   ├─ polyvec_basemul_acc_montgomery (compute secret)
   ├─ poly_invntt_tomont (convert back)
   ├─ poly_sub (subtract from ciphertext)
   ├─ poly_tomsg (decode message)
   ├─ hash_G (re-derive encryption coins)
   ├─ indcpa_enc (re-encrypt to verify)
   ├─ verify (check if ciphertext matches)
   └─ cmov (select correct key)
```

## Files to Create/Modify

| File | Action | Content |
|------|--------|---------|
| **cbd.c / cbd.h** | Create | CBD noise sampling functions |
| **verify.c / verify.h** | Create | Constant-time compare/move |
| **symmetric-shake.c / symmetric.h** | Create | Keccak/SHAKE wrappers |
| **fips202.c / fips202.h** | Create | Keccak/SHA3 math |
| **indcpa.c / indcpa.h** | Create | IND-CPA encryption matrix ops |
| **kem.c / kem.h / api.h** | Create | KEM top-level wrapper |
| **ntt.c / ntt.h** | Update | Add invntt(), basemul(), fqmul() |
| **poly.c / poly.h** | Update | Add all poly ops |
| **polyvec.c / polyvec.h** | Update | Add invntt, basemul, compress, add, reduce |
| **project.c** | Update | Call all functions sequentially |
| **Makefile** | Update | Add new source files to build |

## Important Numbers

- **KYBER_K = 2** (Kyber-512 security level)
- **q = 3329** (the modulus, a prime number)
- **n = 256** (polynomial degree - each polynomial has 256 coefficients)
- **eta1 = 3** (noise parameter for Kyber-512)
- **eta2 = 2** (always 2)
