# Masked KyberNTT Design

Date: 2026-07-31

## Goal

Create a **masked** version of the existing KyberNTT leakage-simulation example so that
the power traces no longer leak the secret coefficient vector `s`. The masked version
uses additive (2-share) masking exploiting the linearity of the Kyber NTT:

    NTT(s) = NTT(s_0) + NTT(s_1)   where   s = s_0 + s_1 (mod q)

The two versions must live in **different directories** and both remain runnable so
their behavior can be compared:

- **Unmasked** — the untouched Kyber example (existing): `elmo/projects/Examples/KyberNTT/`
- **Masked** — new: `elmo/projects/Examples/KyberNTTMasked/`

The masked and unmasked versions are expected to be **different in the trace domain**
(the masked one must not correlate with the secret), while their **computed NTT
results** (printed data) must be identical, because both compute the same correct
`NTT(s)`.

## Security goal

The true secret `s` must never be processed inside the recorded trace window
(`starttrigger()` ... `endtrigger()`). Only the random-looking shares `s_0` and `s_1`
are handled in the window.

## Files / layout

New directory `elmo/projects/Examples/KyberNTTMasked/`, copied from `KyberNTT/`.

Changed files:
- `project.c` — modified for masking. The original unmasked logic is kept in the file as
  comments (so the repo's version stays visible/auditable inside the new file).
- `projectclass.py` — new class `KyberNTTMaskedSimulation`, same input format as the
  unmasked one (it still writes the secret `s`; the device generates its own random share).

Unchanged (copied verbatim): `ntt.c`, `poly.c`, `polyvec.c`, `reduce.c`, `ntt.h`,
`poly.h`, `polyvec.h`, `reduce.h`, `params.h`, `elmoasmfunctionsdef.h`,
`elmoasmfunctionsdef-extension.h`, `elmoasmfunctions.s`, `elmoasmfunctions.o`,
`vector.o`, `project.ld`, `Makefile`.

## C program (`project.c`)

Same structure as the original `project.c`, with the split and double-NTT.

```c
read2bytes(&nb_challenges);
for each challenge:

    // Load the private vector s  (same as unmasked)
    read2bytes(&skpv.vec[j].coeffs[k]);            // for all j, k

    // --- MASKING (all before starttrigger, so s never enters the trace) ---
    for each coefficient (j, k):
        rand2bytes((uint16_t*) &s0.vec[j].coeffs[k]);          // random share
        s0.vec[j].coeffs[k] = barrett_reduce(s0.vec[j].coeffs[k]); // keep in [0, q]
        s1.vec[j].coeffs[k] = skpv.vec[j].coeffs[k] - s0.vec[j].coeffs[k];
        if (s1.vec[j].coeffs[k] < 0) s1.vec[j].coeffs[k] += KYBER_Q; // s1 in [0, q]

    starttrigger();
    // Leaking operations now only involve the shares:
    polyvec_ntt(&s0);
    polyvec_ntt(&s1);
    // Recombine: NTT(s) = NTT(s_0) + NTT(s_1) (mod q)
    for each coefficient (j, k):
        skpv.vec[j].coeffs[k] = barrett_reduce(
            (int16_t)(s0.vec[j].coeffs[k] + s1.vec[j].coeffs[k]));
    endtrigger();

    // Print the (recombined) results, same as unmasked, outside the trace window
    print2bytes(&skpv.vec[j].coeffs[k]);                       // for all j, k

endprogram();
```

### Why the split is done mod q (critical detail)

If `s_1 = s - s_0` is allowed to wrap in 16-bit arithmetic, the recombination is off by
`m * NTT(2^16 mod q)` because `2^16 ≢ 1 (mod 3329)`. Keeping both shares in `[0, q]`
guarantees `s_0 + s_1 ≡ s (mod q)` with no int16 overflow, so the recombined result is
exactly `NTT(s)`.

`poly_ntt` already Barrett-reduces its output (`poly.c:16-20`), so each share's NTT output
is in `[0, q]` and the final single add + `barrett_reduce` is safe.

`#include "reduce.h"` must be added to `project.c` for `barrett_reduce`.

## Python wrapper (`projectclass.py`)

- Class name: `KyberNTTMaskedSimulation` (so `get_simulation('KyberNTTMaskedSimulation')` works).
- `KYBER_K = 2`, `KYBER_N = 256` (same as unmasked).
- `get_binary_path()` returns `'project.bin'`.
- `set_input_for_each_challenge()` writes the secret `s` exactly like the unmasked class.
- `get_test_challenges()` and `get_random_challenges()` copied from the unmasked class.

The only change from the unmasked `projectclass.py` is the class name (and this
documented intent).

## Build

Run `make` inside `elmo/projects/Examples/KyberNTTMasked/` (already done once via
`arm-none-eabi-gcc` 14.3). `project.bin` is produced for ELMO.

## Verification

Both checks run from a Python script using the same set of challenge secrets.

1. **Correctness** — `get_printed_data()` for masked vs unmasked must be **identical**
   (both reduce mod q to the same range, both compute `NTT(s)`).

2. **Defense (traces differ, secret does not leak)** — for the same challenges, Pearson
   correlation of every trace column against `HW(s)` (a standard CPA-style model over
   all `s` coefficients):
   - Unmasked: strong correlation peaks at the leakage points.
   - Masked: correlation ≈ 0 everywhere.
   - Additionally, the raw masked and unmasked traces are not the same (they are
     expected to differ — that is the defense).

Deliverable: a `verify_masking.py` script in the repo root that runs both checks and
prints/saves the correlation plots.

## Out of scope

- More than 2 shares.
- Changing the unmasked `KyberNTT/` directory in any way.
- Making the random share uniform-perfectly (barrett_reduce of `rand2bytes` is uniform
  enough for a leakage demo; real Kyber samples shares from a PRF).
