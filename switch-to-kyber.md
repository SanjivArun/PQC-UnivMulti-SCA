# Plan: Switch tvla_analysis.py from KyberNTT to Kyber

## Current state
- **Commit:** `bb1ad7b` - "docs: add BDV21 A2B wiring plan"
- **Branch:** `master` — pushed to `origin/master`
- **Untracked:** `elmo/projects/Examples/KyberMasked/test_masking` only

## Goal
Change `tvla_analysis.py` to use `elmo/projects/Examples/Kyber/` and `elmo/projects/Examples/KyberMasked/` instead of `elmo/projects/Examples/KyberNTT/` and `elmo/projects/Examples/KyberNTTMasked/`.

## Class name conflict
Both Kyber and KyberMasked define `KyberKEMSimulation`. Need to rename the masked one to `KyberKEMMaskedSimulation`.

## Changes

### Step 1 — Rename class in `KyberMasked/projectclass.py`
- Rename `KyberKEMSimulation` → `KyberKEMMaskedSimulation`

### Step 2 — Update `tvla_analysis.py` paths (lines 55-63)
- `UNMASKED_BIN`: `KyberNTT` → `Kyber`
- `MASKED_BIN`: `KyberNTTMasked` → `KyberMasked`
- VERSIONS list: `KyberNTTSimulation` → `KyberKEMSimulation`, `KyberNTTMaskedSimulation` → `KyberKEMMaskedSimulation`

### Step 3 — Update docstring/header (lines 1-34)
- `KyberNTT` → `Kyber` (and `KyberNTTMasked` → `KyberKEMMasked` where appropriate)

### Step 4 — Update argparse description, plot titles and remaining KyberNTT references
- Line 554: argparse description `"KyberNTT"` → something generic
- Lines 378-385: `fig.suptitle` references to `"KyberNTT"` → `"Kyber"`
