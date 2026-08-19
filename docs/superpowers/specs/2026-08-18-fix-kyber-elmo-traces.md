# Fix Kyber/KyberMasked ELMO Traces

Date: 2026-08-18

## Problem

ELMO hangs indefinitely when running `Kyber/project.bin` and `KyberMasked/project.bin`.

### Root Cause

The compiled `main()` in `Kyber/project.bin` contains an **infinite loop** at address `0x080017b0`:

```asm
80017b0: e7fe  b.n  80017b0  ; infinite loop
```

The compiler optimized away almost all of the `project.c` body. The disassembly shows:
1. Read `nb_challenges` from input
2. If `nb_challenges == 0`, jump to `endprogram`
3. If `nb_challenges > 0`: set up `kgbuf`, call `sha3_512`, call `indcpa_keypair_derand`
4. **Infinite loop** — the rest of the KEM (encapsulation, decapsulation, printing) is stripped

This happens because:
- `project.c` prepares test data from the loop counter `num_challenge` (never read from input)
- `print2bytes()` writes to a memory-mapped register (`0xE0000000`) that the compiler doesn't recognize as having side effects
- The compiler determines the loop body produces no externally visible output and optimizes it to a single iteration
- After the single iteration, the remaining code (enc/dec/prints) is unreachable dead code, producing the infinite loop

**Contrast with KyberNTT**: `KyberNTT/project.c` reads per-challenge NTT coefficients via `read2bytes()` — the compiler must generate the read because it's a side effect. This prevents the optimizer from stripping the loop body.

## Solution

Rewrite `project.c` for both `Kyber/` and `KyberMasked/` to follow the KyberNTT pattern:
1. Read `nb_challenges` from input
2. For each challenge: read per-challenge secret data via `read2bytes()`
3. Wrap the target operation in `starttrigger()/endtrigger()`
4. Print results via `print2bytes()` after `endtrigger()`

This ensures the compiler cannot optimize away the loop body because `read2bytes()` and `print2bytes()` are recognized as having side effects.

## Files to Change

### `elmo/projects/Examples/Kyber/project.c`
- Replace with a version that reads a secret key from input and runs `crypto_kem_dec`
- Follow KyberNTT structure: read challenge, trigger, operate, print, endprogram
- Keep all existing includes (`kem.h`, `indcpa.h`, `symmetric.h`)

### `elmo/projects/Examples/KyberMasked/project.c`
- Same structure as Kyber version
- Use `masked_indcpa_dec` instead of `indcpa_dec` for the decapsulation
- Include `masked_indcpa.h` and `masked_poly.h`

### `elmo/projects/Examples/KyberMasked/projectclass.py`
- No changes needed — already uses `KyberKEMMaskedSimulation` with correct challenge format

### `elmo/projects/Examples/Kyber/projectclass.py`
- No changes needed — already uses `KyberKEMSimulation` with correct challenge format

## Challenge Format

Each challenge provides a Kyber secret key (`KYBER_INDCPA_SECRETKEYBYTES` bytes = `KYBER_POLYVECBYTES` = 768 bytes at K=2).

The `set_input_for_each_challenge` method already writes `(KYBER_K, KYBER_N)` coefficients — 512 `uint16_t` values = 1024 bytes. The secret key is 768 bytes. The mapping: read 512 coefficients, pack them into the secret key format expected by `polyvec_frombytes`.

Actually, the existing `set_input_for_each_challenge` writes raw polynomial coefficients. The `crypto_kem_dec` function expects a packed secret key. The project.c should read the coefficients and use `polyvec_frombytes`-style unpacking, or read the packed form directly.

**Simplest approach**: Read `KYBER_SECRETKEYBYTES` bytes per challenge (1952 bytes at K=2), same format as the existing `set_input_for_each_challenge` writes coefficients. The `project.c` will read coefficients and pack them into the secret key structure.

Actually, the cleanest approach mirrors KyberNTT exactly: read `(K, N)` coefficients and pass them directly to the operation. For Kyber KEM, the operation is `crypto_kem_dec` which needs a full secret key. We'll read the secret key bytes directly.

## Build & Test

1. `make -C elmo/projects/Examples/Kyber clean && make`
2. `make -C elmo/projects/Examples/KyberMasked clean && make`
3. Run ELMO on each binary with 1 challenge, verify traces are produced
4. Run `tvla_analysis.py` with small trace count to verify end-to-end

## What Stays the Same

- All masking code in `masked_poly.c`, `masked_a2b.c`, `masked_indcpa.c` — untouched
- Makefiles — unchanged
- `projectclass.py` files — unchanged (already use correct simulation class names)
- `tvla_analysis.py` — already points to `Kyber/` and `KyberMasked/`
