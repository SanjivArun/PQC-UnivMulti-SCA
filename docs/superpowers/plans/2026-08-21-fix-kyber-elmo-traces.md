# Fix Kyber ELMO Trace Generation (produce 256/group traces)

Date: 2026-08-21
Status: EXECUTING

## Goal
Make `tvla_analysis.py` (default) produce **512 traces** (256 fixed + 256
random) for **both** `KyberKEMSimulation` (unmasked) and
`KyberKEMMaskedSimulation` (masked), so the TVLA table prints
ground-truth / univariate / multivariate leakage counts.

Symptom (user): **ELMO hangs / infinite loop** when running
`elmo/projects/Examples/Kyber/`, so 0 traces are produced.

## Diagnosis

### What "does not compile" is NOT
- The C build is **green**: `make` exits 0, produces `project.bin`
  (36 KB). Only 2 **non-fatal** `memcpy.c` warnings; they are fatal only
  under `-Werror` (not set). RAM budget is fine (64 KB; `.bss`=0,
  `.data`=256 B; big arrays are ~4.7 KB on the stack).

### What it ACTUALLY is — two faces of one root cause
1. **Hang / infinite loop (user's symptom, primary).** The compiler
   optimizes away `project.c`'s loop body -> a `b.n` self-loop
   (documented in `2026-08-18-fix-kyber-elmo-traces.md`, addr
   `0x080017b0`). ELMO never terminates -> 0 traces.
2. **Segfault (seen on the current `-O0` binary, secondary).** Baseline
   repro crashed in the ELMO tool's own code:
   `SEGV @0x68 in gettracelength (fixedvsrandom.h:72) <- write32
   (elmo.c:570) <- execute (3456) <- run (4124)`. Same class of bug; the
   loop body is not "externally visible" to ELMO.

### Root cause (unifying)
The leaky work is derived from `num_challenge` and constants, **not from
the input stream**, and `print2bytes()` writes to an MMIO register
(`0xE0000000`) the compiler/tool don't track as a side effect. So the
compiler can remove the loop body (-> hang) or leave a garbage
half-optimized body (-> tool segfault). This is exactly why **KyberNTT
works** (it reads per-challenge coefficients via `read2bytes()` and so
keeps the loop).

## Experiment plan (controlled bisection)

### Phase 0 - Reproduce (read-only)
- Run `./elmo <abs>/project.bin` (plain + `elmo-asan2`).
- Run `tvla_analysis.py --nb-fixed 4 --nb-random 4 --verbose`.
- Record: **hang / segfault / N traces**. Baseline to compare fixes.

### Control FIRST - KyberNTT
- Run the **known-good** `KyberNTT` binary (1.2 KB).
  - **Produces traces ->** tool/toolchain/`.ld`/RAM are healthy; the
    fault is **Kyber-specific** (the code pattern). Proceed Phase 1.
  - **Also fails ->** fault is the **prebuilt `elmo` / environment**, not
    Kyber code. Jump to Fix #6; do not bisect code.

### Phase 1 - Minimal `project.c` (the "comment out main" step)
- Replace body with the template:
  `read2bytes(&nb)` -> `for: starttrigger(); endtrigger();` ->
  `endprogram()`, nothing else; `make`; run `elmo`.
  - Produces N traces -> toolchain OK; bisect forward.
  - Still fails -> tool/environment problem.

### Phase 2 - Re-add one block at a time
Starting from the working minimal version, after each step
`make` + `elmo`, the first step that re-triggers the failure is the
culprit:
1. read secret via `readbyte` loop (current `sk_inner`)
2. keygen: `hash_g` + `indcpa_keypair_derand`
3. `crypto_kem_enc_derand`
4. `starttrigger()` / `crypto_kem_dec` / `endtrigger()`
5. the `print2bytes` loops

### Phase 3 - Fix
Apply the fix matching Phase 2's culprit (below).

### Phase 4 - Verify 256
- Full `tvla_analysis.py` (256 fixed + 256 random = 512) for **both**
  `KyberKEMSimulation` and `KyberKEMMaskedSimulation`.
- Confirm the comparison table prints for both.

## All the possible fixes (in priority order)

1. **Make each loop iteration side-effect-visible (root cause).** Drive
   the leaky op from **input-read data, not `num_challenge`** - the
   **KyberNTT pattern**. Ensure `project.o` uses `-O0 -fno-inline`
   (it does); ensure the leaky op consumes the **read-back** bytes.
2. **Realign challenge <-> read sizes.** `projectclass.py` writes
   `K*N=512` `uint16` (1024 B)/challenge, but `project.c` reads
   `KYBER_INDCPA_SECRETKEYBYTES` bytes via `readbyte` - a mismatch can
   feed `gettracelength` garbage. Match the two.
3. **`print2bytes` / MMIO `0xE0000000`.** Confirm prints happen strictly
   **after** `endtrigger()` (current code already does); test removing
   prints to isolate.
4. **RAM / `.ld`.** Big arrays are ~4.7 KB on the stack (fits 64 KB);
   verify, or move to a controlled section if a version grows.
5. **Build hygiene (only if a stricter toolchain is the real
   "does not compile").** Fix/drop `memcpy.c` warnings; fix Makefile
   `$(HEADER)` -> `$(HEADERS)`; make `clean` also remove `vector.o` and
   `elmoasmfunctions.o`.
6. **ELMO tool itself (only if the KyberNTT control fails).** The
   prebuilt `elmo` may be the problem - repair/replace it; do not chase
   Kyber code.

## Files in scope
- `elmo/projects/Examples/Kyber/project.c` (primary)
- `elmo/projects/Examples/Kyber/Makefile` (Phase 5 hygiene)
- `elmo/projects/Examples/KyberMasked/project.c` (mirror the fix)
- `elmo/projects/Examples/Kyber/projectclass.py` (only if sizes realigned)
- **Unchanged:** masking sources, `tvla_analysis.py`, `.ld`.

## Success criterion
`tvla_analysis.py` (default) yields **512 traces** for both versions and
prints the ground-truth / univariate / multivariate table.

## Results so far

### Phase 0 (reproduce) + control
- **Kyber** via real pipeline (`get_simulation('KyberKEMSimulation').run()`):
  prints "running Kyber with 10 challenges...", then **HANGS**, 0 traces
  produced. Matches user symptom + `2026-08-18` doc (compiler strips loop
  body -> `b.n` self-loop).
- **KyberNTT control** (`kybernttsimrun.py`, 10 random challenges):
  **produces 10 traces**, shape `(10, 84187)`. -> tool/toolchain/`.ld`/RAM
  are healthy; fault is **Kyber-specific**.
- ASan repro on Kyber (stale input): emulation ran, then SEGV in
  `gettracelength` `feof(NULL)` on a missing `trace00001.trc` — a tool
  NULL-check bug that surfaces once trace generation fails; not the root
  cause, but it masks the true hang when input mismatches.
- Both `Kyber` and `KyberNTT` project classes write the **same** 512
  `uint16` (1024 B) challenge; `Kyber/project.c` reads a different byte
  count via `readbyte` (size mismatch = prime suspect for the hang).

## Bisection results (2026-08-22)

- Build is GREEN (only non-fatal `memcpy.c` warnings). "does not compile" =
  a misnomer; the real failure is **trace generation**.
- **elmo CAN produce traces**: the known-good **KyberNTT** `project.c`
  dropped into the Kyber project produced **10 traces in ~32s** through the
  real pipeline (fixed-vs-random found 668 leaky points). => the tool,
  toolchain, `.ld`, and 64K RAM are all healthy.
- Minimal read (`512 read2bytes` into a `polyvec`, empty trigger) = **10
  traces, fast**. => reading the challenge is NOT the problem.
- `readbyte(50)` works; the full original hangs.
- Removing the **prints** => still hangs. Removing **`crypto_kem_dec`**
  => still hangs. => the hang is **upstream** of `crypto_kem_dec`.
- `static` buffers (.bss 4160) + `-O2` (text 35840 -> 11424) => **still
  hangs** via the reliable pipeline. => not a main-stack overflow; the hang
  is inside the **`crypto_kem_*` emulation** (keypair / SHAKE / X-wing
  path), i.e. a function that is slow-or-infinite in the emulator.
- Caveat: the last crypto bisection steps were **confounded** — the shared
  `elmo/elmo-tool/input.txt` is written by whichever `get_simulation`
  resolves the (duplicate) simulation name, and several runs wrote 0 bytes,
  so elmo took **0 challenges** and exited (not a hang). The **single
  reliable** result is via the pipeline `kyber_run.py` (correct input).

## Conclusion
- **Answer to "can ELMO produce 256 traces?": YES** (elmo works; 10 traces
  proven; ~256 would be ~850s at the KyberNTT pattern's rate).
- The **full-Kyber `project.c`** (keypair+enc+crypto_em_dec) **hangs
  reliably** even with static buffers and `-O2`; the cause is inside the
  emulated crypto (most likely the SHAKE256/X-wing keypair path, which the
  working polyvec_ntt path does not exercise).
- `Makefile` currently changed to `-O2`; `project.c` restored to the
  original full version (with `static` buffers) but it hangs.
