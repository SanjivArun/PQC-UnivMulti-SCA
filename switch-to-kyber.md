# Plan: Switch tvla_analysis.py from KyberNTT to Kyber

## Current state
- **Commit:** `06cfb2a` - "Restructure project.c with 3 trigger blocks (KeyGen/Encap/Decap) for Kyber and KyberMasked"
- **Branch:** `master` — pushed to `origin/master`
- **Untracked:** `elmo/projects/Examples/KyberMasked/test_masking`

## Goal
Change `tvla_analysis.py` to use `elmo/projects/Examples/Kyber/` and `elmo/projects/Examples/KyberMasked/` instead of `elmo/projects/Examples/KyberNTT/` and `elmo/projects/Examples/KyberNTTMasked/`.

## Class name conflict
Both Kyber and KyberMasked define `KyberKEMSimulation`. Renamed the masked one to `KyberKEMMaskedSimulation`.

## Completed: Steps 1–4 — Switch tvla_analysis.py from KyberNTT to Kyber
All 4 original steps are complete in **`55251e6`**:
- Step 1 — Renamed `KyberKEMSimulation` → `KyberKEMMaskedSimulation` in `KyberMasked/projectclass.py`
- Step 2 — Updated `tvla_analysis.py` paths and VERSIONS list from KyberNTT to Kyber/KyberMasked
- Step 3 — Updated docstring/header references
- Step 4 — Updated argparse description and plot titles

## Completed: Steps 5–6 — Fix ELMO compatibility for full Kyber
- Step 5 (`ebc6ef0`) — Increased linker script RAM from 8K → 64K (ROM 64K → 256K) in both `Kyber/project.ld` and `KyberMasked/project.ld` to accommodate full Kyber memory needs
- Step 6 (`06cfb2a`) — Restructured `project.c` in `Kyber/` and `KyberMasked/`:
  - Added `memcpy.c` (minimal implementation with ELMO attributes) to `Kyber/Makefile`
  - Added local buffers (`kgbuf`, `enc_coins`, `enc_kr`) instead of relying on function-scoped vars after loop
  - Split single trigger region into 3 separate trigger blocks: **KeyGen**, **Encap**, **Decap** with matching `endtrigger`
  - Added output printing for public key, secret key, shared secret, and ciphertext
  - Restructured `KyberMasked/project.c` identically (without memcpy.c since it uses different code path)

## Remaining
- Verify TVLA runs end-to-end with the new project structure
- Any further ELMO compatibility issues specific to ELMO framework
- Any masking correctness validation for the new trigger structure
