# Full-Kyber ELMO Design

Date: 2026-08-04

## Goal

Extend the existing NTT-only leakage-simulation example (`elmo/projects/Examples/KyberNTT/`)
so that a single `project.c` runs an entire Kyber-512 KEM round-trip:
`crypto_kem_keypair` → `crypto_kem_enc` → `crypto_kem_dec`, and also exercises each of the
individual leaking sub-operations (forward/inverse NTT, pointwise-Montgomery,
compress/decompress, frommsg/tomsg, poly_add/reduce/tomont, constant-time verify/cmov).

Each leaking operation is wrapped in its own `starttrigger()/endtrigger()` pair, so ELMO
produces **one power-trace file per operation per challenge**. A compile-time flag
`KYBER_TRACE_MODE` selects between per-operation windows (`MULTI_WINDOW`, default) and a
single whole-algorithm window (`SINGLE_WINDOW`).

Everything lives inside `elmo/projects/Examples/KyberNTT/`. No new top-level directories.

## Non-goals / scope

- Unmasked only. No masked sibling is built at this phase.
- No forging of leakage. The real-data TVLA runs must NOT use the `--inject-*` flags; those
  remain a detector self-validation mechanism only.
- K is fixed at 2 (Kyber-512) to fit the 8K RAM region. Building at K>=3 fails fast.

## Backward compatibility

- `KyberNTTSimulation` and the existing `project.bin` (NTT-only) keep working unchanged, so
  `verify_masking.py`, `tvla_analysis.py`, `kybernttsimrun.py`, `test_tvla.py` continue to run.
- A new class `KyberFullSimulation` (added to the same `projectclass.py`) drives a new binary
  `project-full.bin` (built from the same `project.c` compiled with `-DKYBER_FULL`).

## Files / layout (all under `elmo/projects/Examples/KyberNTT/`)

### Vendored from `pq-crystals/kyber` `ref` branch (round3)

Headers (new): `api.h`, `indcpa.h`, `kem.h`, `cbd.h`, `verify.h`, `symmetric.h`,
`fips202.h`, `randombytes.h`.

Headers (replace with ref supersets): `params.h`, `poly.h`, `polyvec.h`, `ntt.h`.
The replaced headers expose MORE functions while still providing the existing NTT-only API
used by `KyberNTTSimulation`.

Sources (new): `cbd.c`, `verify.c`, `indcpa.c`, `kem.c`, `symmetric-shake.c`, `fips202.c`.

Sources (replace with ref supersets): `poly.c`, `polyvec.c`, `ntt.c`. `reduce.c` stays.

### Local files

- `randombytes.c` — ELMO-backed `randombytes()` using `readbyte()` (see below).
- `project.c` — dual `main()`, `#ifdef KYBER_FULL`.
- `projectclass.py` — adds `KyberFullSimulation`.
- `Makefile` — adds `project-full.bin` target (compiles with `-DKYBER_FULL`).

### Namespacing

`KYBER_K` stays 2. `KYBER_NAMESPACE(s)` = `pqcrystals_kyber512_ref_##s` and
`FIPS202_NAMESPACE(s)` = `pqcrystals_kyber_fips202_ref_##s`. `KYBER_90S` is not defined, so
only the SHAKE path (`symmetric-shake.c` + `fips202.c`) compiles. Call sites in `project.c`
use the unprefixed `crypto_kem_keypair/enc/dec` names, which macro-resolve to the prefixed
names — the standard KyberRef idiom.

## C program (`project.c`)

Two `main()` functions chosen by `#ifdef KYBER_FULL`:

- `KYBER_FULL` undefined → NTT-only `main()` (the current `KyberNTT/project.c` body,
  unchanged) → `project.bin`.
- `KYBER_FULL` defined → full-Kyber `main()` → `project-full.bin`.

`KYBER_TRACE_MODE` (`MULTI_WINDOW` default / `SINGLE_WINDOW`) selects the window layout of
the full-Kyber `main()`.

### `randombytes` contract

`randombytes(uint8_t *out, size_t outlen)` fills `out` from the ELMO input stream:

```c
void randombytes(uint8_t *out, size_t outlen) {
  for (size_t i = 0; i < outlen; i++) readbyte(&out[i]);
}
```

Bytes consumed per challenge (drives the challenge format):
- keypair: 64 B (`d||z`)
- enc: 32 B (`m`)
- dec: 0 B (falls back to `z` already stored in `sk` during keypair)
- total = **96 B per challenge**.

### Full-Kyber `main()` per challenge (multi-window)

Per-challenge seed material (96 B) is read once via `randombytes`. Big buffers are `static`
globals to fit 8K RAM at K=2.

N = 14 leaky windows + 1 empty trailer = **15 traces per challenge** (multi-window mode).

| # | Operation (in window) | Printed after `endtrigger` (32 B) |
|---|-----------------------|-----------------------------------|
| 1 | `crypto_kem_keypair(pk, sk)` | `pk[0..31]` |
| 2 | `crypto_kem_enc(ct, ss_enc, pk)` | `ct[0..31]` |
| 3 | `crypto_kem_dec(ss_dec, ct, sk)` | `ss_dec` (full) |
| 4 | `polyvec_ntt(&s_a)` | first 16 coeffs |
| 5 | `polyvec_invntt_tomont(&s_a)` | first 16 coeffs |
| 6 | `polyvec_basemul_acc_montgomery(&s_p, &s_a, &s_b)` | first 16 coeffs |
| 7 | `polyvec_compress(s_buf, &s_a)` | first 32 B compressed |
| 8 | `polyvec_decompress(&s_a, s_buf)` | first 16 coeffs |
| 9 | `poly_compress(s_buf, &s_p)` | first 32 B compressed |
| 10 | `poly_decompress(&s_p, s_buf)` | first 16 coeffs |
| 11 | `poly_tomsg(s_msg, &s_p)` | `s_msg` (full) |
| 12 | `poly_frommsg(&s_p, s_msg)` | first 16 coeffs |
| 13 | `poly_add` + `poly_reduce` + `poly_tomont` | first 16 coeffs |
| 14 | `verify(...)` + `cmov(...)` | verdict + 31 zero pad |
| 15 | empty trailer | (nothing printed inside window) |

After the trailer's `endtrigger()`, print the full dump outside any trigger flag:
`pk` (800) + `sk` (1632) + `ct` (768) + `ss_enc` (32) + `ss_dec` (32) = **3264 B**.

Printing inside the trailer (or anywhere) happens BEFORE its own `endtrigger()` only for
sub-op labels; the full dump is printed strictly AFTER `endtrigger()`, so it contributes no
power samples (the trigger flag gates power sampling; `printbyte` only writes to a file —
see `elmo/elmo-tool/elmo.c` MMIO handler for `0xE0000000`).

Per-challenge printed bytes = `14*32 + 3264 = 3712 B`.

### Correctness assertion

After window 3 (decaps): `assert(memcmp(ss_enc, ss_dec, 32) == 0)`.

### SINGLE_WINDOW mode

Windows 1..14 collapse into one `starttrigger()/endtrigger()`; the trailer stays empty.
Produces 2 traces per challenge (roll-up + trailer). Different binary (`KYBER_TRACE_MODE`
flag) — used for a coarse whole-algorithm leakage profile.

## Printed-data segmentation

`get_printed_data(per_trace=True)` divides the stream equally by `nb_traces`, which does NOT
match the uneven per-window layout. So `KyberFullSimulation` implements a custom segmenter:

```python
def get_printed_data_per_window(self, challenge_index):
    # flat stream via get_printed_data(per_trace=False)
    # per-challenge block = 14*32 labels + 3264 dump
    # returns [{label_i (32B)}, ..., {dump (3264B)}] for the challenge
```

Real bytes, derived from the actual binary output — no injection.

## RAM budget (K=2)

Statics: `pk`+`sk`+`ct`+`ss_enc`+`ss_dec` (3264 B) + `s_a` (1024) + `s_b` (1024) + `s_p` (512)
+ `s_msg` (32) + `s_buf` (up to 640, reused) ≈ 6.5 KB. Stack peak inside `indcpa_*` ~3 KB.
Total < 8 KB with margin. Build-time guard: `#if KYBER_K != 2 #error`.

## Real-data TVLA protocol (to document, not implemented in this phase)

1. Group A = N fixed-seed challenges; Group B = N random-seed challenges.
2. Run both via `KyberFullSimulation.run()` — NO `--inject-*` flags.
3. Apply `tvla.py` Welch per-point and Hotelling per-window at order 1 and 2.
4. Report whatever is found. "Multivariate strictly wins" is an empirical result, honestly
   attributed to the real traces, possibly zero for the unmasked binary.

## Verification of completion

- `make -C elmo/projects/Examples/KyberNTT/ clean` then `make` builds BOTH `project.bin` and
  `project-full.bin` cleanly.
- Old behavior preserved: `verify_masking.py` still passes.
- New behavior: `KyberFullSimulation` with 3 random challenges yields
  `(3*15, L)` traces and `3` per-window segmenters each of length `14*32 + 3264`.
- Correctness assertion holds for all challenges (no decaps failure).