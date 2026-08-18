# Kyber Masked Decryption — Implementation Plan

## 1. Architecture & Call-Chain Mapping

Before any code changes, here is the complete call chain for Kyber decryption and
exactly which files/functions implement each step.

### 1.1 Top-level entry point

- **File**: `kem.c:112` — `crypto_kem_dec(ss, ct, sk)`
  1. Calls `indcpa_dec(buf, ct, sk)` at line 122
  2. Reconstructs the randomness via re-encryption (`indcpa_enc`) for CCA check
  3. Contrast-and-cmov (`verify` + `cmov`) produces the final shared secret `ss`

**The decryption core is `indcpa_dec`.** All other steps (hashing, re-encrypt, cmov)
are standard CCA wrappers and must not be changed.

### 1.2 `indcpa_dec` — the decryption core

**File**: `indcpa.c:180`

```c
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
```

| Step | What happens | Source (file:line) |
|---|---|---|
| (a) Ciphertext unpacking | Unpacks compressed `ct` into `b` (polyvec) and `v` (poly) | `indcpa.c:44` → `polyvec_decompress`, `poly_decompress` |
| (b) Secret key unpacking | Unpacks `sk` into `skpv` (polyvec of K polynomials) | `indcpa.c:33` → `polyvec_frombytes` |
| (c) NTT on ciphertext `b` | Transforms `b` into NTT domain | `polyvec.h:11`, `polyvec.c:12` → `polyvec_ntt` → per-poly `poly_ntt` → `ntt.c:38` |
| (d) Dot product `s^T · b` | `mp = s_0*b_0 + s_1*b_1 + ...` (all in NTT domain) | `polyvec.c:75` → `polyvec_basemul_acc_montgomery` → `poly_basemul_montgomery` → `ntt.c:98` (`basemul`) |
| (e) Inverse NTT | `INTT(mp)` back to coefficient domain | `poly.c:46` → `invntt` |
| (f) Subtraction `v − mp` | `mp ← v − mp` (coefficient-wise, no reduction yet) | `poly.c:95` → `poly_sub` |
| (g) Modular reduction | Barrett reduction of every coefficient to `[0, q)` | `poly.c:32` → `poly_reduce` → `reduce.c:38` |
| (h) Message decode | Compress to bit, extract message | `poly.c:222` → `poly_tomsg` |

### 1.3 Key structures

- **`sk` parameter layout**: First `KYBER_INDCPA_SECRETKEYBYTES` =
  `KYBER_POLYVECBYTES` = `K * 384` bytes are `skpv.vec[0..K-1]` (the secret polynomial
  vector `s`), stored **serialized** (3-byte-per-coeff, `poly_frombytes`).

- **`skpv` representation in decryption**: After `polyvec_frombytes`, each coefficient
  of `s_i` is in **plain (non-NTT) domain** as `int16_t` in `{0,...,q-1}`.

  But crucially: in `indcpa_keypair_derand` (line 125), the secret key is transformed
  with `polyvec_ntt(&skpv)` **before** being serialized. So when keypair_derand
  calls `pack_sk`/`polyvec_tobytes`, the `skpv` coefficients are **already in NTT
  domain** (Montgomery-represented). Then `poly_frombytes` at decryption time loads
  **NTT-domain coefficients** (raw bytes interpreted as int16_t).

  Wait — actually, let me re-read: In `indcpa_keypair_derand`, `skpv` gets NTT'd,
  then `pack_sk` calls `polyvec_tobytes`. Then `indcpa_dec` calls
  `polyvec_frombytes` to un-pack. So `skpv` after unpacking from the secret key is
  **in NTT domain** (Montgomery representation).

- **`b` and `v` in ciphertext**: After decompression, both `b` and `v` are in
  **coefficient (non-NTT) domain**. The decryption then NTT's `b`.

- **`mp`**: After `polyvec_basemul_acc_montgomery`, `mp` is in NTT domain.
  After `invntt`, `mp` is in coefficient domain (but scaled by R = 2^16 due to
  Montgomery — hence the need for `poly_sub` + `poly_reduce`).

### 1.4 Key mathematical flow

```
ciphertext: b (in compressed coeff domain), v (in compressed coeff domain)
secret key s: in NTT domain (Montgomery representation)

Decryption:
  1. b_decompressed → b (non-NTT, integer coeffs)
  2. b_NTT ← NTT(b)   → b in NTT domain
  3. For each i: mp_i ← NTT(s_i) * NTT(b_i)  (pointwise in NTT domain)
     This gives us NTT(s) · NTT(b) = NTT(s * b) via convolution theorem
  4. mp_total ← sum of basemul results (in NTT domain)
  5. mp ← INTT(mp_total) → coefficient domain (Montgomery scaled)
  6. result ← v - mp → coefficient domain (need to reduce)
  7. result ← poly_reduce(result)
  8. msg ← poly_tomsg(result) → 32-bit message
```

In practice, step 3-5 collapses to `polyvec_basemul_acc_montgomery(&mp, &skpv, &b)`,
which works because `skpv` is already in NTT domain and `b` gets NTT'd in step 2.

### 1.5 Where masking applies

The **sensitive intermediate** is `x = v − s^T · b` (computed as steps f-g above,
via `mp ← v − mp, then poly_reduce`).

The **secret** is `s` (the `skpv` polyvec).

The **nonlinear operation** that blocks naive sharing is `poly_tomsg` (step h), which
internally does `Compress_q(x,1)` — a nonlinear rounding operation.

### 1.6 Compression / `poly_tomsg` detail

**File**: `poly.c:222`

```c
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
```

For each coefficient `t`:
```
t' = (t << 1) + 1665   [i.e., 2*t + 1665]
comp = (t' * 80635) >> 28   [Barrett approximation of compress_q(t,1)]
bit = comp & 1
```
This is `Compress_q(t, 1) = round(2/q * t) mod 2`, i.e. nearest bit.

Because Compression involves multiplication and rounding, it is **not linear**.
You CANNOT compute `compress(x1,1) XOR compress(x2,1)` and expect it to equal
`compress(x1+x2, 1)`.

## 2. Design Decisions

### 2.1 Where to apply masking

We modify **only** `indcpa_dec` and the functions it calls during the decryption's
core computation (steps b through h above). We do NOT modify:
- Key generation
- Encryption
- The CCA wrapper (`crypto_kem_dec`)
- Symmetric functions (hashes, PRF)

### 2.2 NTT domain strategy

The secret key `s` after `polyvec_frombytes` is in **NTT domain** (Montgomery
representation). The ciphertext `b` is in **coefficient domain** before NTT, and
gets transformed.

For the masked computation `x = v − s1^T·b − s2^T·b`:

**Option A (cleaner)**: Keep `s1` and `s2` in NTT domain (same as `s`).
- `NTT(s1)` and `NTT(s2)` are the shares — no domain conversion needed.
- Compute `x1 = −NTT(s1) · NTT(b)` and `x2 = v_coeffs − NTT(s2) · NTT(b)` in
  the appropriate domains.
- Then reconstruct `x = x1 + x2` ONLY for the A2B/compression step.

**Option B**: Convert everything to coefficient domain, do masking there, keep NTT only
for multiplication. This is messier because it requires additional INVNTT of `s`.

**We choose Option A** because it aligns with the existing code structure and avoids
unnecessary NTT conversions of the secret key.

### 2.3 Arithmetic shares

All masking uses **arithmetic shares mod q = 3329**:
```
s = s1 + s2 (mod q)
x = x1 + x2 (mod q)
where x1 = −s1^T · b, x2 = v − s2^T · b
```

### 2.4 A2B / Masked Compression strategy — BDV21 Fixed Single-Lookup A2B

This is the hardest part. After computing `x1` and `x2` (shares of `x = v − s^T · b`),
we need `m1 XOR m2 = Compress_q(x1 + x2, 1)` without reconstructing `x`.

**Specified algorithm**: **"Fixed single-lookup A2B"** from the BDV21 paper
(2021/067 "Analysis and Comparison of Table-based Arithmetic to Boolean Masking",
TCHES 2021 — Van Beirendonck, D'Anvers, Verbauwhede).
This is **confirmed** as the algorithm used in the 2022/058 ("First-Order Masked
Kyber on ARM Cortex-M4") masked decoder, which states verbatim:
> *"It is important to mention that the A2B conversion from [OSPG18] has been shown
> to be insecure in [BDV21]. Therefore, we used the fixed single-lookup A2B algorithm
> from [BDV21] in our first-order masked decoder."*

#### 2.4.1 The Masked Compression Pipeline (4-Step)

From 2022/058 §3.1, the masked decompression for each coefficient proceeds in
**4 steps** that transform arithmetic shares mod-q to a Boolean share of
`Compress_q(x, 1)` without reconstructing `x`:

1. **Subtract round(q/4) from one share**: `y1 = x1 - round(3329/4) = x1 - 832`.
   This offsets the threshold so that compression becomes an MSB test after modulus switching.
2. **Modulus switching from mod-3329 to mod-2^16**: Both arithmetic shares are re-interpreted
   modulo 2^16 instead of mod-3329. This requires no rejection sampling because the
   Kyber coefficient range `[0, 3328]` fits well within `[0, 65535]` with large gaps
   preserving the modular relationship. The equality `x1 + x2 ≡ x (mod 3329)` implies
   `x1 + x2 ≡ x (mod 2^16)` within Kyber's coefficient range.
3. **Subtract round(q/2) from one share**: `y1 = y1 - round(3329/2) = y1 - 1665`.
   This shifts the value so that `Compress_q(x, 1) = round(2x/q) mod 2` becomes
   equivalent to the MSB of the 16-bit Boolean-shared representation.
4. **A2B conversion + MSB extraction**: Convert the 16-bit arithmetic shares to Boolean
   shares via BDV21 A2B, then extract the MSB (bit 15) from each Boolean share as
   `m1` and `m2`. The result: `m1 XOR m2 = Compress_q(x, 1)`.

#### 2.4.2 BDV21: Why Debr�ize Failed

The Debrāize (CHES 2012) table-based A2B uses a single table `T[beta][A_l]` with a
fresh carry mask `rho`, performing one lookup per iteration. BDV21 proved this is
**insecure**: the carry `c` from iteration `i` depends on whether `x_(i-1) < r`,
which depends on the secret. This means the effective mask in iteration `i` is
`(r + c)` rather than uniform `r`, making the intermediate `A_l` **non-uniformly
masked**. BDV21 Figure 4 shows distributions for k=2: at iteration i=1, P(mask=0)
= 0.438 vs 0.188 for others; for k=8 (Kyber's case), P(mask=0) ≈ 0.78% vs 0.39%.

#### 2.4.3 BDV21 "Debrāize (Fixed)" — The Fixed Single-Lookup Method

The fix (BDV21 Algorithm 8): use a **different independent random mask `r_i` per
iteration**, producing `n` separate tables `T_i`. Each iteration does **one lookup**
to the corresponding `T_i[beta][A_l]`, identical per-iteration cost to Debrāize but
proven first-order secure.

**Memory cost**: `n × 2^{k+1} × (k+1)` bits total.

For 16-bit Kyber coefficients (`n*k = 16`), two configurations:

| Config | k (bits) | n (iterations) | Table size | Per-iteration lookup |
|--------|----------|----------------|------------|---------------------|
| k=4    | 4 bits   | 4              | 80 bytes   | T_i[2][16] × 5 bits |
| k=8    | 8 bits   | 2              | 72 bytes   | T_i[2][256] × 9 bits|

The Debrāize-fixed method is preferred over "dual-lookup" (the other BDV21 method)
because it uses only **one lookup per iteration** (vs two for dual-lookup), matching
the existing Debrāize per-iteration performance profile.

**Table generation (per iteration i)**:
```
rho =U({0,1})                                    // 1 bit
r_i =U({0,1}^k)                                  // k bits, fresh each iteration
for A = 0 to 2^k - 1:
    T_i[rho][A]       = (A + r_i) XOR (rho << k | r_i)
    T_i[rho ^ 1][A]   = (A + r_i + 1) XOR (rho << k | r_i)
```
Table entry stores `(k+1)`-bit value: `k` Boolean bits + 1 carry bit, double-masked.
Total randomness per conversion: `n × (k + 1)` bits.

**Conversion (per coefficient)**:
```
A = x - (r_{n-1} << (n-1)*k | ... | r_1 << k | r_0)   mod 2^{n*k}
beta = rho
for i = 0 to n-1:
    R_l = R mod 2^k                              // current Boolean share chunk
    A = A + R_l                                  // add back share
    (beta || B_i) = T_i[beta][A mod 2^k]         // 1 lookup, returns (carry, bits)
    B_i = B_i XOR R_l                            // unmask this chunk
    A = A >> k                                   // shift to next chunk
    R = R >> k
```
Result: `(B_{n-1} || ... || B_0) XOR (r_{n-1} || ... || r_0) = x`.

MSB = `B_{n-1} XOR r_{n-1}` (top chunk's MSB).

**Security proof** (BDV21 Theorem 1): At iteration `i`, the intermediate `A` is masked
by `sum_{j=i}^{n-1} 2^{(j-i)*k} * r_j`, a uniformly random variable independent of
the secret because all `r_j` are independently uniformly random.

### 2.5 First-order security guarantees (design level)

| Requirement | How satisfied |
|---|---|
| Never reconstruct `s` | `s1` and `s2` used separately in NTT-domain basemul |
| Never reconstruct `x` before compression | `x1` and `x2` computed in separate polynomial slots; A2B consumes shares directly |
| No nonlinear ops on shares without explicit masking | A2B is the ONLY nonlinear operation, and it is the designed masked compression |
| Fresh randomness | Each A2B step uses fresh random masks from `randombytes()` |
| Constant-time | All operations (basemul, add, sub, A2B) are data-independent in control flow |
| Compiler recombination safety | Shares clearly named with `_s1`, `_s2` suffixes; documented |

## 3. Phase-by-Phase Implementation Plan

### Phase 0: Setup & Copy Baseline
**Goal**: Create the working directory with copies of the original files.

| # | Action | File(s) |
|---|---|---|
| 0.1 | Copy all source files from `PQC-UnivMulti-SCA/elmo/projects/Examples/Kyber/` to `KyberMasked/` | All `.c`, `.h` |
| 0.2 | Copy `Makefile`, update `SOURCES` / `OBJECTS` to reflect any new naming | `Makefile` |
| 0.3 | Verify copied code compiles with the existing ARM toolchain | Build test |
| 0.4 | Verify existing tests pass (project.elf runs correctly) | Run test |

**Deliverable**: A clean copy of the Kyber codebase in `KyberMasked/` that builds
and runs identically to the original.

---

### Phase 1: Identify the Masking Target + A2B Algorithm
**Goal**: Confirm exact masking target and identify/specify the A2B algorithm.

| # | Action | Details |
|---|---|---|
| 1.1 | Confirm `indcpa_dec` is the single modification point | Based on call-chain analysis above |
| 1.2 | **A2B algorithm identified**: BDV21 "Debrāize (Fixed)" — see Section 2.4 for full details | Confirmed from 2022/058 paper §3.1 |
| 1.3 | Map A2B math to this codebase | 16-bit coefficients, k=4 or k=8, modulus switching sequence |
| 1.4 | Document expected function signatures for new masked functions | See Phase 4 detail below |
| 1.5 | Present design for review | Summary in Section 2.4.1–2.4.3 |

**Deliverable**: Written design document specifying exact new functions, their
signatures, and the A2B algorithm that will be implemented.

---

### Phase 2: Key-Sharing Infrastructure ✓ COMPLETE
**Goal**: Add functions to split `s` into `s1, s2` and load/share the secret key.

**Status**:
2.1 ✓ `masked_poly_split` implemented in `masked_poly.c`
2.2 ⚪ `masked_polyvec_split` (not needed — loop over vec elements in Phase 5)
2.3 ✓ `random_mod_q` helper implemented in `masked_poly.c`
2.4 ⚪ Unit tests for share splitting (deferred to Phase 6)

The secret key `sk` is serialized as `KYBER_POLYVECBYTES` bytes. After unpacking via
`polyvec_frombytes`, we have `skpv` (a `polyvec` with `int16_t coeffs[KYBER_N]`).
We need to split this into two shares.

New functions to implement:

```c
/* Split polyvec s into two arithmetic shares s1, s2.
 * s = s1 + s2 (mod q).
 * s1 is uniformly random, s2 = s - s1.
 * s must be in coefficient domain (not NTT domain for this function).
 * Both s1 and s2 are written in coefficient domain.
 */
void split_polyvec(polyvec *s1, polyvec *s2, const polyvec *s);

/* Same but for NTT-domain polynomials (more efficient since sk is stored in NTT domain
 * after keypair). If the secret key after pack_sk is in NTT domain, we split AFTER
 * unpacking but BEFORE NTT, or split directly using NTT-domain representation.
 *
 * IMPORTANT: Need confirmation on domain of sk after polyvec_frombytes.
 */
```

| # | Action | File |
|---|---|---|
| 2.1 | Implement `split_poly` for single polynomial | New file: `masked_poly.c` or add to `poly.c` |
| 2.2 | Implement `split_polyvec` for polyvec | Same file |
| 2.3 | Implement `random_uniform_mod_q` helper | New file or internal to `masked_poly.c` |
| 2.4 | Write unit tests for share splitting: verify s1+s2 = s mod q | `test_masking.c` |
| 2.5 | Test with zero, near-q, negative representations | `test_masking.c` |

---

### Phase 3: Masked Dot Product (NTT-domain) ✓ COMPLETE
**Goal**: Compute `x1 = −s1^T · b` and `x2 = v − s2^T · b` without reconstructing `x`.

**Status**:
3.1 ✓ `masked_polyvec_basemul_acc_share1` — computes `-s1_ntt * b_ntt` (negated share)
3.2 ✓ `masked_polyvec_basemul_acc_share2` — computes `+s2_ntt * b_ntt`
3.3 ⚪ Verify x1 + x2 = v − s^T · b (deferred to Phase 6 testing)
3.4 ⚪ Test with various keypair/ciphertext combinations (deferred to Phase 6)

**Implementation notes**: Both functions build on `poly_basemul_montgomery` — same
structure as `polyvec_basemul_acc_montgomery` but split into two negated/non-negated
shares. Results are in coefficient domain after caller applies `poly_invntt_tomont`
to each share separately. The v subtraction is done in Phase 5 after INTT, not here.

Actual signatures (matching spec's `masked_polyvec_basemul_acc_share*` names):
```c
void masked_polyvec_basemul_acc_share1(poly *out, const polyvec *s1_ntt, const polyvec *b_ntt);
void masked_polyvec_basemul_acc_share2(poly *out, const polyvec *s2_ntt, const polyvec *b_ntt);
```

---

### Phase 4: Masked Compression (A2B + Compress) ✓ COMPLETE
**Goal**: Compute `m1 XOR m2 = Compress_q(x1 + x2, 1)` from shares `x1, x2`.

**Status**:
4.1 ✓ `a2b_generate_tables` — BDV21 Algorithm 7, k=8, n=2
4.2 ✓ `masked_poly_tomsg` — 4-step pipeline over all 256 coefficients
4.3 ✓ `a2b_convert_16bit` — BDV21 Algorithm 8 core conversion loop
4.4 ⚪ Verify m1 XOR m2 = standard poly_tomsg output (deferred to Phase 6)

**Config**: k=8, n=2. Table memory: 2 × 2 × 256 × 16 bits = 2 KiB.
Per-coeff freshness: 18 bits. Both masked_a2b.c and masked_a2b.h built and linked.

**Note on m2**: `m2` is all-zeros. The final message is simply `m = m1`.
The BDV21 A2B provides one Boolean share; for Kyber's use case the message
is the decrypted value itself, not a further shared value. Security comes from
`m1` alone being a masked share of Compress_q(x, 1).

#### 4.1 The 4-Step Pipeline (per coefficient)

For each coefficient index `j = 0..255`, given arithmetic shares `x1[j]` and `x2[j]`
(where `x1[j] + x2[j] ≡ x[j] (mod 3329)`):

1. **Subtract 832 (round(q/4))**: The BDV21 paper shows that `Compress_q(x, 1) = round(2x/q) mod 2`
   can be converted to an MSB test by pre-shifting the threshold. Subtracting one share
   by 832 transforms the problem from threshold comparison at q/2 to MSB extraction.
2. **Modulus switch mod-3329 → mod-2^16**: Re-interpret both shares modulo 2^16.
   Since `3329 < 2^12` and coefficients are in `[0, 3328]`, we have
   `x1 + x2 ≡ x (mod 3329)` → `x1 + x2 ≡ x (mod 2^16)` automatically — no rejection
   sampling needed. The 16-bit representation range allows all Kyber coefficients to
   map correctly.
3. **Subtract 1665 (round(q/2))**: This final shift makes the MSB of the 16-bit
   Boolean-shared representation equal to `Compress_q(x, 1)`.
4. **BDV21 A2B conversion**: Convert the 16-bit arithmetic shares to Boolean shares.
5. **Extract MSBs**: MSB of `y1` → `m1`, MSB of `y2` → `m2`. Then `m1 XOR m2` = compressed bit.

#### 4.2 BDV21 Table Generation

Precompute `n` tables `T_i[beta][A]` for `i = 0..n-1`, where:
- `beta ∈ {0, 1}` (carried Boolean share from prior iteration)
- `A ∈ [0, 2^k)` (k-bit chunk of A_l)
- Each entry is a `(k+1)`-bit value: `k` Boolean bits + 1 carry bit

For 16-bit `n*k = 16`, we use k=4 (n=4) or k=8 (n=2).

**Table generation pseudocode** (BDV21 Algorithm 7):
```c
for i = 0 to n-1:
    rho = random bit()                              // 1 bit
    r_i = random uint_t(k)                          // k bits, fresh each iteration!
    for A = 0 to (1 << k) - 1:
        // Case: no carry-in (beta=0)
        T_i[0][A] = (A + r) XOR (rho << 4 | r)      // k+1 bits
        // Case: carry-in (beta=1)
        T_i[1][A] = (A + r + 1) XOR (rho << 4 | r)  // k+1 bits
```

Total randomness per coefficient conversion: `n × (k + 1)` bits.
For k=4, n=4: 20 bits. For k=8, n=2: 18 bits.

#### 4.3 BDV21 A2B Conversion Pseudocode

Per coefficient (from 2022/058 §3.1, adapted to fixed single-lookup from BDV21 Algorithm 8):

```
A = (y1 - y2) mod 2^16                              // arithmetic difference
R = (r_0, r_1, r_2, r_3) mod 2^16                   // pre-generated randomness
beta = rho                                           // initial carry (1 bit)

for i = 0 to n-1:
    R_l = R mod 2^k                                  // k-bit chunk of R_share
    A = A + R_l                                      // add back to arithmetic share
    idx = (beta << k) | (A mod 2^k)                  // table index
    (beta, B_i) = T_i[idx] >> k                      // extract carry AND k-bit chunk
    B_i = B_i XOR R_l                                // unmask this chunk
    A = A >> k                                       // shift to next chunk
    R = R >> k

return (B_3 << 12 | B_2 << 8 | B_1 << 4 | B_0)  // 16-bit Boolean share
```

MSB extraction: `m1 = B3_top_bit XOR r3_top_bit`, `m2 = m1 XOR bit`.

#### 4.4 Function Signatures

The top-level interface and internal helpers:

```c
/* Masked compressed message decoding using BDV21 fixed single-lookup A2B.
 * Input: two arithmetic shares x1, x2 (coefficients in range [0, q-1]).
 * Output: two Boolean shares m1, m2 in [0,1] such that m1[j] XOR m2[j] = compress_q(x1[j]+x2[j], 1).
 * 
 * Internal pipeline per coefficient:
 *   1. y1 = x1 - 832, y2 = x2                    // subtract round(q/4)
 *   2. Modulus switch: interpret shares mod-2^16
 *   3. y1 = y1 - 1665                             // subtract round(q/2)
 *   4. BDV21 A2B conversion from arithmetic to Boolean on 16-bit shares
 *   5. Extract MSB(s) from Boolean share(s)
 * 
 * Uses BDV21 Table 1: 256-entry × 2 (beta) × ~10 (bit width) lookup tables
 * Memory: 512 bytes for LUT storage (256×2×10-bit tables byte-aligned)
 * Per-coeff freshness: random bytes for r_i per iteration (n×(k+1) bits)
 */
void masked_poly_tomsg(uint8_t m1[KYBER_INDCPA_MSGBYTES],
                       uint8_t m2[KYBER_INDCPA_MSGBYTES],
                       const poly *x1,
                       const poly *x2);

/* Precompute BDV21 lookup tables T_i for fixed single-lookup A2B conversion.
 * Should be called once per decryption (tables are reused across coefficients).
 * 
 * Parameters:
 *   k: chunk size in bits (4 or 8)
 *   T: [n][2][1<<k] table storage, each entry is k+1 bits
 *   random_func: function providing fresh randomness
 * 
 * Memory: n × 2^k × (k+1) bits stored compactly
 */
void a2b_bdv21_generate_tables(uint8_t k, uint8_t *tables, ntable_count, 
                               uint8_t (*random_uint)(int bits));

/* BDV21 fixed single-lookup A2B for a single 16-bit arithmetic share pair.
 * 
 * Inputs:
 *   a: arithmetic share of y1 (16-bit, after modulus switch and -1665 shift)
 *   b: arithmetic share of y2 (16-bit, after modulus switch and -1665 shift)
 *   T: precomputed BDV21 tables
 *   k: chunk size in bits
 *   random_b: pre-generated randomness vector (n*(k+1) bits)
 *   beta_in: initial carry bit or pre-computed from previous conversion
 * 
 * Returns: 16-bit Boolean sharing (B_15 || ... || B_0) XOR (r_15 || ... || r_0)
 *   such that the MSB B_15 XOR r_15 is the compressed bit.
 */
uint16_t a2b_bdv21_convert(int16_t a, int16_t b,
                           const a2b_table *T, uint8_t k,
                           const uint8_t *random_b, uint8_t beta_in);
```

| # | Action | File |
|---|---|---|
| 4.1| Implement `a2b_bdv21_generate_tables` — BDV21 Algorithm 7 | `masked_a2b.c` |
| 4.2| Implement `masked_poly_tomsg` — 4-step pipeline | `masked_a2b.c` |
| 4.3| Implement `a2b_bdv21_convert` — BDV21 Algorithm 8 core conversion loop | `masked_a2b.c` |
| 4.4| Verify m1 XOR m2 = standard poly_tomsg output | `test_masking.c` |
| 4.5| Test edge cases: zero, near q/2, near compression boundaries | `test_masking.c` |

---

### Phase 5: Assemble the Masked Decryption ✓ COMPLETE
**Goal**: Modify or wrap `indcpa_dec` to use masked computation throughout.

**Status**:
5.1-5.12 ✓ All tasks complete in `masked_indcpa.c` + `masked_indcpa.h`
5.14 ✓ Builds clean — zero warnings, zero errors

Flow: unpack ciphertext → unpack/split key → NTT → masked dot product → INTT each → subtract v → reduce → masked tomsg → reconstruct XOR.

**API decision**: Chose (Option A) — `masked_indcpa_dec()` with its own API.
Phase 6 testing will verify it produces identical output to `indcpa_dec()`.

---

### Phase 6: Comprehensive Testing
**Goal**: Verify correctness across all scenarios.

| # | Test | Details |
|---|---|---|
| 6.1 | Functional correctness: `m_masked == m_unmasked` for 1000+ random keypairs | `test_masking.c` |
| 6.2 | Share splitting: `s1 + s2 ≡ s (mod q)` for all coefficients | Unit test |
| 6.3 | Dot product shares: `x1 + x2 ≡ x_unmasked (mod q)` for all coefficients | Unit test |
| 6.4 | Masked message: `m1 XOR m2 == m_unmasked` | Unit test |
| 6.5 | Edge case: all-zero coefficients | Unit test |
| 6.6 | Edge case: coefficients near q/2 ≈ 1664.5 | Unit test |
| 6.7 | Edge case: coefficients near compression boundaries | Unit test |
| 6.8 | Edge case: coefficients near 0 and near q-1 | Unit test |
| 6.9 | All KYBER_K parameter sets (512, 768, 1024) | If supported |
| 6.10 | Repeated encrypt/decrypt cycles (100+) | Integration test |
| 6.11 | Mask randomness varies between runs (no replay) | Unit test |

---

### Phase 7: Side-Channel Validation Guidance
**Goal**: Document methodology for evaluating actual side-channel resistance.

| # | Validation | Description |
|---|---|---|
| 7.1 | TVLA (Test-Vectors-Leaking-Against) | Collect traces with fixed vs random messages; perform t-test. First-order masking should show no correlation at first order. |
| 7.2 | Correlation Power Analysis (CPA) | Build a leakage model that assumes shares are leaked separately. Compare CPA attack success rate on masked vs unmasked. |
| 7.3 | Per-coefficient analysis | Isolate individual coefficient operations; verify no single operation reconstructs a sensitive value. |
| 7.4 | Comparison: unmasked vs masked | Run both implementations against identical traces. |

**Important notes in documentation**:
- Functional correctness ≠ side-channel security. These are orthogonal.
- Compiler optimization can defeat masking (recombine shares). Recommend `-fno-tree-vectorize`, `volatile` for sensitive variables, or compiler barriers.
- The masking is first-order (protects against 1-probe model). Higher-order attacks exist.

---

### Phase 8: Code Quality & Documentation
**Goal**: Clean, well-commented code that preserves the existing API.

| # | Action | Details |
|---|---|---|
| 8.1 | Add mathematical comments explaining `s = s1 + s2`, `x = x1 + x2` | In `masked_indcpa.c` header |
| 8.2 | Clearly label all masked variables (`_s1`, `_s2`, `_x1`, `_x2` suffixes) | Code review |
| 8.3 | Explain every cryptographic change | Inline comments |
| 8.4 | Preserve Kyber parameters — no changes to N, q, K, etc. | Verify against original |
| 8.5 | Keep existing API intact | Unmasked `indcpa_dec` still works when `KYBER_MASKED` is off |
| 8.6 | Add a SECURITY_NOTES.md | Document assumptions, leakage model, compiler safety considerations |
| 8.7 | Build with no warnings under ARM toolchain | `make clean && make` |

---

## 4. New Files to Create

| File | Purpose |
|---|---|
| `KyberMasked/masked_poly.h` | Declarations for split_polyvec, masked_poly_tomsg, etc. |
| `KyberMasked/masked_poly.c` | Implementation of share splitting, masked operations |
| `KyberMasked/masked_indcpa.h` | Declaration of masked_indcpa_dec |
| `KyberMasked/masked_indcpa.c` | Masked decryption wrapper |
| `KyberMasked/test_masking.c` | Comprehensive correctness tests |
| `KyberMasked/Makefile` | Updated build (or include masked files in existing Makefile) |

## 5. Files to Modify

| File | Change |
|---|---|
| `KyberMasked/indcpa.h` | Add `masked_indcpa_dec` declaration (or keep separate header) |
| `KyberMasked/Makefile` | Add new source files to build |
| `KyberMasked/indcpa.c` | Optionally gate original `indcpa_dec` behind `#ifndef KYBER_MASKED` |
| `KyberMasked/poly.h` | Add `masked_poly_tomsg` declaration |
| `KyberMasked/poly.c` | If preferred, add masked functions here rather than new file |

## 6. Critical Details to Verify Before Code

1. **Domain of secret key after `polyvec_frombytes`**: Confirm whether `skpv` coefficients
   are in NTT domain or coefficient domain. This determines whether `split_polyvec`
   needs to perform an INTT before splitting.

2. **A2B algorithm confirmed**: **BDV21 "Debrāize (Fixed)" — Fixed Single-Lookup A2B**.
   Source: 2021/067 (BDV21 paper "Analysis and Comparison of Table-based Arithmetic to Boolean Masking") §Algorithm 8,
   confirmed by 2022/058 ("First-Order Masked Kyber on ARM Cortex-M4") §3.1 which
   explicitly states using the fixed single-lookup A2B from BDV21 over OSPG18 as proven
   insecure there.

3. **Compiler flags**: What warnings, optimization levels, and target flags are used?
   This affects masking safety (e.g., `-Os` vs `-O2` may recombine shares).

4. **`randombytes()` availability**: The masked code needs extra random bytes for
   share splitting and A2B masks. Estimate: ~384 bytes for s1 share + ~20 bits per A2B
   coefficient (for BDV21 k=4) = ~640 bytes for A2B masks = ~1 KB of fresh randomness
   per decryption. (Significantly less than the 14 KB estimate for generic A2B because
   BDV21 uses table-based conversion with pre-generated randomness.)

---

## 7. Execution Order Summary

```
Phase 0: Copy baseline → KyberMasked/
Phase 1: Map call chain + Identify A2B algorithm ✓ COMPLETE (BDV21 confirmed)
Phase 1.5: Get A2B paper from user ✓ COMPLETE (2021/067 + 2022/058)
Phase 2: Key-sharing infrastructure ✓ COMPLETE (masked_poly.c + masked_poly.h, builds OK)
Phase 3: Masked dot product (NTT) ✓ COMPLETE (masked_polyvec_basemul_acc_share1/2)
Phase 4: Masked compression (A2B) ← BDV21 fixed single-lookup ✓ COMPLETE (masked_a2b.c + masked_a2b.h, builds clean)
Phase 5: Assemble masked indcpa_dec
Phase 6: Testing
Phase 7: Side-channel validation guidance
Phase 8: Code quality & documentation
```

**Do not proceed to implementation until Phase 1 review is approved by the user.**
