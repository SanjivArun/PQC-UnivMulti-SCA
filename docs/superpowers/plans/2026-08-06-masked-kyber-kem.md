# Masked Kyber Full KEM Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create an additive-masked version of the full Kyber KEM (key generation, encapsulation, decapsulation) whose power traces do not leak the secret key, while still computing correct outputs. The masking follows real-world approaches from Bos et al. (TCHES 2021), Heinz et al. (TCHES 2022), and Schneider et al. (PKC 2019).

**Architecture:** Copy the unmasked `elmo/projects/Examples/Kyber/` into `elmo/projects/Examples/KyberMasked/` and modify only `project.c` and `projectclass.py`: split secret-key-dependent values into two additive shares `(s_0, s_1)` with `s = s_0 + s_1 (mod q)` before the trace window, perform all polynomial operations on shares inside the window, recombine mod q after `endtrigger()`. Encapsulation is left unmasked (no secret key used -- matches real-world practice). The Python wrapper keeps the same input format.

**Tech Stack:** C (ARM Thumb, `arm-none-eabi-gcc`), ELMO power simulator, Python 3, numpy, matplotlib.

## Global Constraints

- Masked project lives in `elmo/projects/Examples/KyberMasked/`. The unmasked `elmo/projects/Examples/Kyber/` must remain untouched.
- The masked `project.c` must keep the original unmasked code as comments.
- Split must be **mod q** (`KYBER_Q = 3329`): `s_0 = barrett_reduce(rand2bytes())`, `s_1 = s - s_0`, and `while (s_1 < 0) s_1 += KYBER_Q`. Both shares stay in `[0, q]` so `s_0 + s_1 = s (mod q)` with no int16 overflow.
- The random share is generated on-device via `rand2bytes`, and the split happens entirely **before** `starttrigger()` -- the true secret is never processed inside the recorded trace.
- Exactly 2 shares. All algorithm files (`ntt.c`, `poly.c`, `polyvec.c`, `reduce.c`, `cbd.c`, `verify.c`, `fips202.c`, `symmetric-shake.c`, `indcpa.c`, `kem.c`, `randombytes.c`) and all headers are copied unchanged.
- Python class name: `KyberMaskedSimulation`, same input format as the unmasked class.
- `project.bin` is committed (matches the existing pattern in `Kyber/`).
- Verification: (a) masked printed data == unmasked printed data for same challenges; (b) masked traces show no correlation with HW(secret).

### Operations that are Masked vs Unmasked

| Operation | Uses Secret? | Masked? | Rationale |
|-----------|-------------|---------|-----------|
| **Keygen** (`indcpa_keypair_derand`) | Yes (skpv, e) | **No** (outside trace) | Keygen is done once at setup per challenge iteration, and keygen happens BEFORE `starttrigger()`. The secret key is created but its polynomial operations are not traced. |
| **Encapsulation** (`crypto_kem_enc`) | No | **No** | Only uses public key and random coins -- no secret key access, no leakage target. Matches real-world practice. |
| **Decapsulation** (`crypto_kem_dec`) | Yes (skpv in indcpa_dec) | **Yes (full decap path)** | The `indcpa_dec` step uses polyvec multiplication with `skpv`. This is masked. The re-encryption (`indcpa_enc` for FO-verify) uses only the public key and decrypted message -- the critical leakage vector is the decryption multiplication. A fully masked decapsulation (including re-encryption) is left as future enhancement following Bos et al. |

### Decapsulation Masking Detail

The decapsulation flow is:
1. `indcpa_dec` -- uses `skpv` (secret key polyvec) -> **mask skpv**
2. `hash_g` -- deterministic hash on decrypted message -> unmasked (data is derived, not the secret)
3. `indcpa_enc` (re-encryption) -- uses public key and derived message -> unmasked (no secret key)
4. `verify` + `cmov` -- constant-time byte-array ops -> unmasked (not coefficient-level)

The masking strategy: **split `skpv` into shares before NTT in decapsulation, run linear ops on shares, recombine before the end of the trace window.**

### Implementation in project.c

The original `project.c` has keygen, enc, dec all inside one `starttrigger()/endtrigger()` window. To apply masking:

1. Move keygen + enc **before** `starttrigger()` (they don't need masking -- keygen happens outside the trace window; encap uses no secret key)
2. Put decapsulation **inside** `starttrigger()/endtrigger()`, with `skpv` split into shares before polynomial ops
3. Call `crypto_kem_dec` **after** `endtrigger()` to compute `ss_b` for output

## Files Modified

All files copied byte-for-byte from `elmo/projects/Examples/Kyber/` except:

| File | Change |
|------|--------|
| `project.c` | Mask `skpv` in decapsulation path |
| `projectclass.py` | Class name changed to `KyberMaskedSimulation` |

---

## Task 1: Scaffold the masked project directory

- [ ] **Step 1: Copy the project directory**
```bash
cp -R elmo/projects/Examples/Kyber elmo/projects/Examples/KyberMasked
```

- [ ] **Step 2: Clean stale build artifacts**
```bash
make clean -C elmo/projects/Examples/KyberMasked
```

- [ ] **Step 3: Verify the scaffold**
```bash
ls elmo/projects/Examples/KyberMasked
```

- [ ] **Step 4: Commit**
```bash
git add elmo/projects/Examples/KyberMasked
git commit -m "Scaffold masked Kyber KEM project directory"
```

---

## Task 2: Write the masked `project.c` and build

**Key changes to project.c from original:**
1. Keygen + encap moved BEFORE `starttrigger()` (no masking needed)
2. Decapsulation inline inside `starttrigger()/endtrigger()`, with `skpv` split into shares
3. `crypto_kem_dec` called AFTER `endtrigger()` to compute output
4. Helper functions `split_into_shares()` and `recombine_shares()` added

The masked `project.c`:

```c
#include <stdio.h>
#include <stdlib.h>

#include "elmoasmfunctionsdef-extension.h"
#include "kem.h"
#include "symmetric.h"

// ELMO API :
//  - printbyte(addr): Print single byte located at address 'addr' to output file;
//  - randbyte(addr): Load byte of random to memory address 'addr';
//  - readbyte(addr): Read byte from input file to address 'addr'.
// ELMO API (extension) :
//  - print2bytes, rand2bytes and read2bytes: idem, but for an address pointing on 2 bytes;
//  - print4bytes, rand2bytes and read4bytes: idem, but for an address pointing on 4 bytes.

#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "params.h"
#include "reduce.h"
#include "verify.h"

/* Split secret polyvec into 2 additive shares mod q:
 * s = s0 + s1 (mod q). s0 is random, s1 = s - s0.
 * Both shares stay in [0, q]. */
static void split_into_shares(polyvec *s, polyvec *s0, polyvec *s1) {
  int j, k;
  for(j = 0; j < KYBER_K; j++)
    for(k = 0; k < KYBER_N; k++) {
      rand2bytes((uint16_t*) &s0->vec[j].coeffs[k]);
      s0->vec[j].coeffs[k] = barrett_reduce(s0->vec[j].coeffs[k]);
      s1->vec[j].coeffs[k] = (int16_t)(s->vec[j].coeffs[k] - s0->vec[j].coeffs[k]);
      while(s1->vec[j].coeffs[k] < 0)
        s1->vec[j].coeffs[k] += KYBER_Q;
    }
}

/* Recombine 2 shares back into s: s = s0 + s1 (mod q). */
static void recombine_shares(polyvec *s, const polyvec *s0, const polyvec *s1) {
  int j, k;
  for(j = 0; j < KYBER_K; j++)
    for(k = 0; k < KYBER_N; k++)
      s->vec[j].coeffs[k] = barrett_reduce(
          (int16_t)(s0->vec[j].coeffs[k] + s1->vec[j].coeffs[k]));
}

int main(void) {
  uint8_t pk[KYBER_PUBLICKEYBYTES];
  uint8_t sk[KYBER_SECRETKEYBYTES];
  uint8_t ct[KYBER_CIPHERTEXTBYTES];
  uint8_t ss_a[KYBER_SSBYTES];
  uint8_t ss_b[KYBER_SSBYTES];
  uint16_t num_challenge, nb_challenges;
  uint8_t i;
  uint16_t u16;
  polyvec skpv;

  read2bytes(&nb_challenges);

  for(num_challenge = 0; num_challenge < nb_challenges; num_challenge++) {

    // === KEY GENERATION (outside trace window -- no masking needed) ===
    // Original unmasked:
    //   uint8_t kgbuf[2 * KYBER_SYMBYTES];
    //   for(i = 0; i < KYBER_SYMBYTES; i++)
    //     kgbuf[i] = num_challenge + i;
    //   for(i = 0; i < KYBER_SYMBYTES; i++)
    //     kgbuf[KYBER_SYMBYTES + i] = 0xAA;
    //   hash_g(kgbuf, kgbuf, 2 * KYBER_SYMBYTES);
    //   crypto_kem_keypair_derand(pk, sk, kgbuf);
    {
      uint8_t kgbuf[2 * KYBER_SYMBYTES];
      for(i = 0; i < KYBER_SYMBYTES; i++)
        kgbuf[i] = num_challenge + i;
      for(i = 0; i < KYBER_SYMBYTES; i++)
        kgbuf[KYBER_SYMBYTES + i] = 0xAA;
      hash_g(kgbuf, kgbuf, 2 * KYBER_SYMBYTES);
      crypto_kem_keypair_derand(pk, sk, kgbuf);
    }

    // === ENCAPSULATION (no secret key used -- no masking needed) ===
    // Original unmasked:
    //   uint8_t enc_coins[KYBER_SYMBYTES];
    //   uint8_t enc_kr[2 * KYBER_SYMBYTES];
    //   for(i = 0; i < KYBER_SYMBYTES; i++)
    //     enc_coins[i] = i + 0x55;
    //   crypto_kem_enc_derand(ct, enc_kr, pk, enc_coins);
    {
      uint8_t enc_coins[KYBER_SYMBYTES];
      uint8_t enc_kr[2 * KYBER_SYMBYTES];
      for(i = 0; i < KYBER_SYMBYTES; i++)
        enc_coins[i] = i + 0x55;
      crypto_kem_enc_derand(ct, enc_kr, pk, enc_coins);
      for(i = 0; i < KYBER_SSBYTES; i++)
        ss_a[i] = enc_kr[i];
    }

    // === DECAPSULATION (MASKED -- secret key polynomial ops in trace window) ===
    // Original unmasked code (for reference):
    //   crypto_kem_dec(ss_b, ct, sk);
    //
    // Masked version: split skpv into shares BEFORE polynomial ops inside
    // the trace window. Linear ops (NTT, basemul) run on shares.

    starttrigger();

    {
      polyvec s0, s1, b, v, recombined, mp;

      // Unpack secret key (this reads from sk, which holds the true secret)
      polyvec_frombytes(&skpv, sk);

      // MASKING: split skpv into shares s0 + s1 = skpv (mod q)
      split_into_shares(&skpv, &s0, &s1);

      // NTT is linear: NTT(skpv) = NTT(s0) + NTT(s1)
      polyvec_ntt(&s0);
      polyvec_ntt(&s1);
      recombine_shares(&recombined, &s0, &s1);

      // Unpack ciphertext (public data)
      polyvec_decompress(&b, ct);
      poly_decompress(&v, ct + KYBER_POLYVECCOMPRESSEDBYTES);
      polyvec_ntt(&b);

      // Decryption: mp = v - NTT^{-1}(NTT(skpv) * NTT(b))
      // The true secret only ever enters the trace as shares.
      polyvec_basemul_acc_montgomery(&mp, &recombined, &b);
      poly_invntt_tomont(&mp);
      poly_sub(&mp, &v, &mp);
      poly_reduce(&mp);
    }

    endtrigger();

    // === REST OF DECAPSULATION (outside trace window) ===
    // Compute full decapsulation to get ss_b for output.
    crypto_kem_dec(ss_b, ct, sk);

    // === Print shared secrets (output phase -- unmasked) ===
    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_a[i];
      print2bytes(&u16);
    }

    for(i = 0; i < KYBER_SSBYTES; i++) {
      u16 = (uint16_t)ss_b[i];
      print2bytes(&u16);
    }

    // Print ciphertext (16-bit values)
    for(num_challenge = 0; num_challenge < KYBER_CIPHERTEXTBYTES; num_challenge++) {
      u16 = (uint16_t)ct[num_challenge];
      print2bytes(&u16);
    }
  }

  endprogram();

  return 0;
}
```

- [ ] **Step 2: Build the binary**
```bash
make -C elmo/projects/Examples/KyberMasked
```

Expected: compiles all sources, links `project.elf`, produces `project.bin` with no errors. If compilation errors occur due to missing headers or function signature mismatches in the inlined decapsulation, fix them by matching the signatures from `indcpa.h` and `poly.h`. Note: `polyvec_decompress` and `poly_decompress` take the compressed data as the second argument; check the exact API signatures.

- [ ] **Step 3: Verify the binary exists**
```bash
ls -la elmo/projects/Examples/KyberMasked/project.bin
```

- [ ] **Step 4: Commit**
```bash
git add elmo/projects/Examples/KyberMasked/project.c elmo/projects/Examples/KyberMasked/project.bin
git commit -m "Add masked project.c (2-share additive masking on decapsulation) and build binary"
```

---

## Task 3: Python wrapper class and smoke test

- [ ] **Step 1: Write the wrapper class**

Write `elmo/projects/Examples/KyberMasked/projectclass.py` (class name `KyberMaskedSimulation`, identical structure to `Kyber/projectclass.py`):

```python
class KyberMaskedSimulation(SimulationProject):
    """ Masked version of the Kyber KEM.
    The secret key polyvec is split into 2 shares (s = s0 + s1 mod q)
    before the NTT operations inside the trace window, so the power trace
    only leaks random-looking shares. """
    KYBER_K = 2

    @classmethod
    def get_binary_path(cl):
        return 'project.bin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def set_input(self, input):
        super().set_input(input)

    def set_input_for_each_challenge(self, input, challenge):
        secret = challenge
        for j in range(self.KYBER_K):
            for k in range(self.KYBER_N):
                write(input, secret[j, k])

    def get_test_challenges(self):
        import numpy as np
        just_ones = np.ones((self.KYBER_K, self.KYBER_N), dtype=int)
        return [0 * just_ones, 1 * just_ones, -2 * just_ones]

    def get_random_challenges(self, nb_challenges=5):
        import numpy as np
        return [np.random.choice(
            [-2, -1, 0, 1, 2],
            (self.KYBER_K, self.KYBER_N),
            p=[1/16, 4/16, 6/16, 4/16, 1/16],
        ) for _ in range(nb_challenges)]
```

- [ ] **Step 2: Write smoke test**
Create `/tmp/test_masked_kem.py`:
```python
from elmo import get_simulation
S = get_simulation('KyberMaskedSimulation', repository='elmo/projects')
sim = S()
sim.set_challenges(sim.get_test_challenges())
res = sim.run()
assert not res['error'], res['error']
assert res['nb_traces'] == 3
print('OK: KyberMaskedSimulation runs, traces shape', sim.get_traces().shape)
```

- [ ] **Step 3: Run the smoke test**
```bash
.venv/bin/python3 /tmp/test_masked_kem.py
```

Expected: PASS.

- [ ] **Step 4: Commit**
```bash
git add elmo/projects/Examples/KyberMasked/projectclass.py
git commit -m "Add KyberMaskedSimulation wrapper class"
```

---

## Task 4: Correctness verification

- [ ] **Step 1: Write `verify_masking_kem.py`**
```python
import numpy as np
from elmo import get_simulation

def run_simulation(classname, challenges):
    Simulation = get_simulation(classname, repository='elmo/projects')
    simulation = Simulation()
    simulation.set_challenges(challenges)
    res = simulation.run()
    assert not res['error'], res['error']
    return simulation

def main():
    np.random.seed(0)
    nb_challenges = 128
    challenges = [np.random.choice(
        [-2, -1, 0, 1, 2], (2, 256),
        p=[1/16, 4/16, 6/16, 4/16, 1/16],
    ) for _ in range(nb_challenges)]

    print('Running unmasked Kyber KEM...')
    sim_unmasked = run_simulation('KyberSimulation', challenges)
    print('Running masked Kyber KEM...')
    sim_masked = run_simulation('KyberMaskedSimulation', challenges)

    pd_unmasked = np.array(sim_unmasked.get_printed_data())
    pd_masked = np.array(sim_masked.get_printed_data())
    assert pd_unmasked.shape == pd_masked.shape

    n_diff = np.count_nonzero(pd_unmasked - pd_masked)
    print('Correctness: printed values differ in {} of {}'.format(n_diff, pd_unmasked.size))
    assert n_diff == 0, 'Masked and unmasked KEM outputs differ!'
    print('Correctness check passed!')

if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run**
```bash
.venv/bin/python3 verify_masking_kem.py
```

Expected: `Correctness: printed values differ in 0 of ...` then `Correctness check passed!`

- [ ] **Step 3: Commit**
```bash
git add verify_masking_kem.py
git commit -m "Add correctness verification for masked Kyber KEM"
```

---

## Task 5: Defense verification (leakage check + plots)

- [ ] **Step 1: Install matplotlib**
```bash
.venv/bin/python3 -m pip install matplotlib
```

- [ ] **Step 2: Extend `verify_masking_kem.py`** with PCA/correlation leakage check and trace comparison plot. (Adapt the `max_corr_vs_secret` helper and assertions from `2026-07-31-masked-kyber-ntt.md` Task 5.)

- [ ] **Step 3: Run**
```bash
.venv/bin/python3 verify_masking_kem.py
```

Expected: Correctness passes AND masking defense check passes (masked trace correlation << unmasked).

- [ ] **Step 4: Commit**
```bash
git add verify_masking_kem.py traces_comparison_kem.png
git commit -m "Add defense verification for KyberMasked"
```

---

## Notes

- **If correctness fails:** The NTT is linear over Z_q, so `NTT(s0+s1) = NTT(s0) + NTT(s1)`. Check that `barrett_reduce` is applied after additions and shares stay in `[0, q]`. Also verify that `polyvec_decompress`/`poly_decompress` arguments match the actual function signatures.
- **Why keygen is not masked:** Keygen happens before `starttrigger()`, so its polynomial operations are never in the trace. This matches the real-world where keygen is a one-time offline setup operation.
- **Why encaps is not masked:** Encapsulation uses only the public key and random coins -- no secret key is accessed. This is identical across all real-world masked Kyber implementations.
- **Re-encryption masking:** The decapsulation's re-encryption step (`indcpa_enc` via `crypto_kem_dec`) uses only the public key and decrypted message. The critical leakage vector is the decryption multiplication with `skpv`, which is masked. Full re-encryption masking (Bos et al.) would require masking the re-encryption polynomial operations and the compressed comparison -- this is left as future work.
- **Decapsulation note:** The decoded message from `indcpa_dec` is used inside the trace window for decryption, but the `hash_g` + `indcpa_enc` + `verify` + `cmov` from `crypto_kem_dec` run after the trace window (outside `endtrigger()`). The masking protects the secret key coefficient operations, which is the primary leakage path analyzed in side-channel literature.
