# Masked KyberNTT Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a 2-share additive-masked version of the KyberNTT leakage-simulation example whose power traces do not leak the secret `s`, while still computing the correct `NTT(s)`.

**Architecture:** Copy the unmasked example into a new directory `elmo/projects/Examples/KyberNTTMasked/` and change only `project.c`: split the secret into two shares `s_0, s_1` (mod q) before the trace window, run `polyvec_ntt` on each share inside the window, recombine mod q, and print the recombined result. The Python wrapper keeps the same input format. Two Python verification checks confirm (a) printed NTT results match the unmasked version exactly and (b) the masked traces show no correlation with `HW(s)`.

**Tech Stack:** C (ARM Thumb, `arm-none-eabi-gcc` 14.3), ELMO power simulator, Python 3.14 (venv `.venv/`), numpy, matplotlib.

## Global Constraints

- Masked project lives in `elmo/projects/Examples/KyberNTTMasked/`. The unmasked `elmo/projects/Examples/KyberNTT/` must remain untouched.
- The masked `project.c` must keep the original unmasked code as comments.
- Split must be **mod q** (`KYBER_Q = 3329`): `s_0 = barrett_reduce(rand2bytes())`, `s_1 = s - s_0`, and `while (s_1 < 0) s_1 += KYBER_Q`. Both shares must stay in `[0, q]` so `s_0 + s_1 ≡ s (mod q)` with no int16 overflow.
- The random share is generated on-device via `rand2bytes`, and the split happens entirely **before** `starttrigger()` — the true secret `s` is never processed inside the recorded trace.
- Exactly 2 shares. `ntt.c`, `poly.c`, `polyvec.c`, `reduce.c` and all headers are copied unchanged.
- Python class name: `KyberNTTMaskedSimulation`, same input format as the unmasked class (writes the secret `s`).
- All Python commands use `.venv/bin/python3` (the venv at repo root, Python 3.14).
- `project.bin` is committed (matches the existing pattern in `KyberNTT/`, whose `.gitignore` only ignores `*.elf`, `*.list`, `*.map`, `*.d`, `*.o`).
- Verification expectations: unmasked `max |corr| vs HW(s)` is high (> 0.4); masked is low (< 0.5 * unmasked and < 0.5).

---

### Task 1: Scaffold the masked project directory

**Files:**
- Create: `elmo/projects/Examples/KyberNTTMasked/` (copy of `elmo/projects/Examples/KyberNTT/`)

**Interfaces:**
- Consumes: nothing.
- Produces: a full copy of the Kyber project that Task 2 modifies.

- [ ] **Step 1: Copy the project directory**

```bash
cp -R elmo/projects/Examples/KyberNTT elmo/projects/Examples/KyberNTTMasked
```

- [ ] **Step 2: Clean stale build artifacts so `project.bin` gets rebuilt from the new C code later**

```bash
make clean -C elmo/projects/Examples/KyberNTTMasked
```

Expected: no errors; removes `.o`, `.elf`, `.bin`, `.list`, `.map`, `.d` files.

- [ ] **Step 3: Verify the scaffold**

```bash
ls elmo/projects/Examples/KyberNTTMasked
```

Expected: contains `project.c`, `projectclass.py`, `Makefile`, `ntt.c`, `poly.c`, `polyvec.c`, `reduce.c`, all `.h` files, `elmoasmfunctions.s`, `elmoasmfunctions.o`, `vector.o`, `project.ld`, `.gitignore`.

- [ ] **Step 4: Commit**

```bash
git add elmo/projects/Examples/KyberNTTMasked
git commit -m "Scaffold masked KyberNTT project directory"
```

---

### Task 2: Write the masked `project.c` and build

**Files:**
- Modify: `elmo/projects/Examples/KyberNTTMasked/project.c`

**Interfaces:**
- Consumes: `elmoasmfunctionsdef-extension.h` (ELMO API: `read2bytes`, `rand2bytes`, `print2bytes`, `starttrigger`, `endtrigger`, `endprogram`), `polyvec.h` (`polyvec`, `polyvec_ntt`), `params.h` (`KYBER_K`, `KYBER_N`, `KYBER_Q`), `reduce.h` (`barrett_reduce`).
- Produces: compiled `project.bin` whose printed data equals the unmasked version's printed data for the same challenges (verified in Task 4).

- [ ] **Step 1: Replace `project.c` with the masked version**

Write the following to `elmo/projects/Examples/KyberNTTMasked/project.c`. The original unmasked body is preserved as comments at the top of the relevant sections.

```c
#include <stdio.h>
#include <stdlib.h>

#include "elmoasmfunctionsdef-extension.h"

// ELMO API :
//  - printbyte(addr): Print single byte located at address 'addr' to output file;
//  - randbyte(addr): Load byte of random to memory address 'addr';
//  - readbyte(addr): Read byte from input file to address 'addr'.
// ELMO API (extension) :
//  - print2bytes, rand2bytes and read2bytes: idem, but for an address pointing on 2 bytes;
//  - print4bytes, rand4bytes and read4bytes: idem, but for an address pointing on 4 bytes.

#include "polyvec.h"
#include "params.h"
#include "reduce.h"

int main(void) {
  uint16_t num_challenge, nb_challenges;
  int j, k;
  polyvec skpv, s0, s1;

  read2bytes(&nb_challenges);
  for(num_challenge=0; num_challenge<nb_challenges; num_challenge++) {

    // Load the private vector s
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++)
        read2bytes((uint16_t*) &skpv.vec[j].coeffs[k]);

    // MASKED VERSION (the original unmasked body is kept below in comments):
    //
    // Original code:
    //   starttrigger(); // To start a new trace
    //   // Do the leaking operations here...
    //   polyvec_ntt(&skpv);
    //   endtrigger(); // To end the current trace
    //
    // Masking: split s into two shares s0, s1 with s = s0 + s1 (mod q).
    // s0 is random, s1 is derived. The split happens BEFORE starttrigger()
    // so that the true secret s never appears inside the recorded trace.
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++) {
        rand2bytes((uint16_t*) &s0.vec[j].coeffs[k]);
        s0.vec[j].coeffs[k] = barrett_reduce(s0.vec[j].coeffs[k]); /* s0 in [0, q] */
        s1.vec[j].coeffs[k] = (int16_t)(skpv.vec[j].coeffs[k] - s0.vec[j].coeffs[k]);
        while(s1.vec[j].coeffs[k] < 0)
          s1.vec[j].coeffs[k] += KYBER_Q; /* s1 in [0, q] */
      }

    starttrigger(); // To start a new trace

    // Leaking operations... now only on the random-looking shares.
    // NTT is linear over Z_q: NTT(s) = NTT(s0) + NTT(s1).
    polyvec_ntt(&s0);
    polyvec_ntt(&s1);

    // Recombine the two shares (mod q)
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++)
        skpv.vec[j].coeffs[k] = barrett_reduce(
            (int16_t)(s0.vec[j].coeffs[k] + s1.vec[j].coeffs[k]));

    endtrigger(); // To end the current trace

    // Print the results of the computation
    for(j=0;j<KYBER_K;j++)
      for(k=0;k<KYBER_N;k++)
        print2bytes((uint16_t*) &skpv.vec[j].coeffs[k]);
  }

  endprogram(); // To indicate to ELMO that the simulation is finished

  return 0;
}
```

- [ ] **Step 2: Build the binary**

```bash
make -C elmo/projects/Examples/KyberNTTMasked
```

Expected: compiles `ntt.o`, `poly.o`, `polyvec.o`, `reduce.o`, `project.o`, links `project.elf`, produces `project.bin` and `project.list`, no errors/warnings.

- [ ] **Step 3: Verify the binary exists**

```bash
ls -la elmo/projects/Examples/KyberNTTMasked/project.bin
```

Expected: a non-empty `project.bin` file.

- [ ] **Step 4: Commit**

```bash
git add elmo/projects/Examples/KyberNTTMasked/project.c elmo/projects/Examples/KyberNTTMasked/project.bin
git commit -m "Add masked project.c (2-share additive masking) and build binary"
```

---

### Task 3: Python wrapper class and smoke test

**Files:**
- Create: `elmo/projects/Examples/KyberNTTMasked/projectclass.py`

**Interfaces:**
- Consumes: `SimulationProject`, `write` (provided by the simulation loader in `elmo/manage.py`), `elmo/projects/Examples/KyberNTTMasked/project.bin`.
- Produces: class `KyberNTTMaskedSimulation(SimulationProject)` — `get_simulation('KyberNTTMaskedSimulation', repository='elmo/projects')` returns it; `run()` + `get_traces()` work.

- [ ] **Step 1: Write the failing test**

Create `/tmp/test_masked_class.py` (outside the repo):

```python
from elmo import get_simulation

S = get_simulation('KyberNTTMaskedSimulation', repository='elmo/projects')
sim = S()
sim.set_challenges(sim.get_test_challenges())
res = sim.run()
assert not res['error'], res['error']
assert res['nb_traces'] == 3
traces = sim.get_traces()
assert traces.shape[0] == 3
print('OK: KyberNTTMaskedSimulation runs, traces shape', traces.shape)
```

- [ ] **Step 2: Run the test to verify it fails**

Run from the repo root:

```bash
.venv/bin/python3 /tmp/test_masked_class.py
```

Expected: FAIL with `SimulationNotFoundError` (no `KyberNTTMaskedSimulation` class yet).

- [ ] **Step 3: Write the wrapper class**

Write the following to `elmo/projects/Examples/KyberNTTMasked/projectclass.py` (identical to the unmasked wrapper except the class name and docstring):

```python
### In this file is defined a Python class to manipulate the simualtion project.
###  - This class must be inherited from th class 'SimulationProject' (no need to import it)
###  - You can use here the function "write(input_file, uint, nb_bits=16)"
###            to write an integer of 'nb_bits' bits in the 'input_file' (no need to import it too).
### To get this simulation class in Python scripts, please use the functions in manage.py as
###  - search_simulations(repository)
###  - get_simulation(repository, classname=None)
###  - get_simulation_via_classname(classname)

class KyberNTTMaskedSimulation(SimulationProject):
    """ Masked version of the KyberNTT example.
    The secret vector s is split into 2 shares (s = s0 + s1 mod q) inside the
    C code, before the trace window, so the power trace only leaks the shares.
    The printed NTT result is the recombined NTT(s). """
    KYBER_K = 2 #k=2 for Kyber512
    KYBER_N = 256 #n=256 for Kyber512

    @classmethod
    def get_binary_path(cl):
        return 'project.bin'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def set_input(self, input):
        """ Write into the 'input' file of ELMO tool
                the parameters and the challenges for the simulation """
        super().set_input(input)

    def set_input_for_each_challenge(self, input, challenge):
        """ Write into the 'input' file of ELMO tool
                the 'challenge' for the simulation """
        secret = challenge

        # Write the secret vector
        for j in range(self.KYBER_K):
            for k in range(self.KYBER_N):
                write(input, secret[j,k])

    def get_test_challenges(self):
        import numpy as np
        just_ones = np.ones((self.KYBER_K, self.KYBER_N), dtype=int)
        return [
             0 * just_ones,
             1 * just_ones,
            -2 * just_ones,
        ]

    def get_random_challenges(self, nb_challenges=5):
        import numpy as np
        return [ np.random.choice(
            [-2, -1, 0, 1, 2],
            (self.KYBER_K, self.KYBER_N),
            p=[1/16, 4/16, 6/16, 4/16, 1/16],
        ) for _ in range(nb_challenges) ]
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
.venv/bin/python3 /tmp/test_masked_class.py
```

Expected: PASS, prints `OK: KyberNTTMaskedSimulation runs, traces shape (3, ...)`.

- [ ] **Step 5: Commit**

```bash
git add elmo/projects/Examples/KyberNTTMasked/projectclass.py
git commit -m "Add KyberNTTMaskedSimulation wrapper class"
```

---

### Task 4: Correctness verification

**Files:**
- Create: `verify_masking.py` (repo root)

**Interfaces:**
- Consumes: `get_simulation` from `elmo`, `KyberNTTSimulation`, `KyberNTTMaskedSimulation`, `project.bin` in both project dirs.
- Produces: script that passes iff masked printed data == unmasked printed data for the same challenges.

- [ ] **Step 1: Write `verify_masking.py` (correctness part)**

Write the following to `verify_masking.py` at the repo root:

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
        [-2, -1, 0, 1, 2],
        (2, 256),
        p=[1/16, 4/16, 6/16, 4/16, 1/16],
    ) for _ in range(nb_challenges)]

    print('Running unmasked KyberNTT...')
    sim_unmasked = run_simulation('KyberNTTSimulation', challenges)
    print('Running masked KyberNTT...')
    sim_masked = run_simulation('KyberNTTMaskedSimulation', challenges)

    # --- Correctness: both versions must compute the same NTT(s) ---
    pd_unmasked = np.array(sim_unmasked.get_printed_data())
    pd_masked = np.array(sim_masked.get_printed_data())
    assert pd_unmasked.shape == pd_masked.shape, (pd_unmasked.shape, pd_masked.shape)
    # barrett_reduce maps zero-residue coefficients to 0 (non-negative inputs)
    # or 3329 (negative multiples of q), so compare modulo q.
    n_diff = np.count_nonzero((pd_unmasked - pd_masked) % 3329)
    print('Correctness: printed values differ in {} of {}'.format(n_diff, pd_unmasked.size))
    assert n_diff == 0, 'Masked and unmasked NTT results differ!'

    print('Correctness check passed!')


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run the correctness check**

```bash
.venv/bin/python3 verify_masking.py
```

Expected: prints both "Running..." lines, then `Correctness: printed values differ in 0 of 131072` (128 challenges × 512 coefficients × 2 bytes) and `Correctness check passed!`. (Takes a few minutes — 128 traces per simulation.)

- [ ] **Step 3: Commit**

```bash
git add verify_masking.py
git commit -m "Add correctness verification for masked KyberNTT"
```

---

### Task 5: Defense verification (leakage check + plots)

**Status: COMPLETE** (see `.git/sdd/task-5-report.md`). Final results with 256 challenges:
correctness 0/131072 coefficients differ mod q; unmasked max |corr| = 1.000; masked
max |corr| = 0.329 (pure noise floor, tracks `~sqrt(2*ln(P)/N)`). `traces_comparison.png` written.

Deviations from the brief below (all necessary):
- **`project.c` leak fix** (also committed with this task): the recombination loop ran
  inside `starttrigger()/endtrigger()`, writing the true NTT(s) into the traced window
  (masked corr 0.961). Moved the recombination to after `endtrigger()` and rebuilt `project.bin`.
- **Correctness byte order**: `print2bytes` emits each uint16 big-endian (high byte first);
  coefficients are rebuilt as `pd[2k]*256 + pd[2k+1]` before the mod-q comparison. The earlier
  "12/131072 differ" finding was this artifact (e.g. bytes `13,1` = 3329 ≡ 0 vs `0,0`).
- **Shared output dir**: both simulations share `elmo/elmo-tool/output/`, so each
  simulation's printed data/traces are captured immediately after its own `run()`.
- `max_corr_vs_secret` skips zero-variance trace columns; plot uses per-subplot x-arrays
  (masked trace 168352 pts vs unmasked 93142 pts); `nb_challenges` raised to 256 so the
  noise floor sits well below the 0.5 threshold.

**Files:**
- Modify: `verify_masking.py` (repo root)
- Modify: `elmo/projects/Examples/KyberNTTMasked/project.c` (recombination moved after `endtrigger()`)
- Modify: `elmo/projects/Examples/KyberNTTMasked/project.bin` (rebuilt)
- Add: `traces_comparison.png`

**Interfaces:**
- Consumes: everything from Task 4, plus `matplotlib` (installed into the venv here).
- Produces: script that prints `Unmasked max |corr| vs HW(s)` and `Masked max |corr| vs HW(s)`, asserts the defense margins, and saves `traces_comparison.png`.

- [ ] **Step 1: Install matplotlib into the venv**

```bash
.venv/bin/python3 -m pip install matplotlib
```

Expected: installs `matplotlib` (and deps) into the venv.

- [ ] **Step 2: Extend `verify_masking.py`**

Add the following helper before `main()`:

```python
def hw(x):
    """Hamming weight of a 16-bit value."""
    x = int(x) & 0xFFFF
    c = 0
    while x:
        c += x & 1
        x >>= 1
    return c


def max_corr_vs_secret(traces, hws):
    """Max |Pearson correlation| between any trace column and HW(secret)."""
    n, npoints = traces.shape
    t = traces - traces.mean(axis=0)
    t_std = t.std(axis=0)
    worst = 0.0
    for c in range(hws.shape[1]):
        h = hws[:, c]
        h_c = h - h.mean()
        denom = n * t_std * h_c.std()
        with np.errstate(invalid='ignore', divide='ignore'):
            corr = (t * h_c[:, None]).sum(axis=0) / denom
        m = np.abs(np.nan_to_num(corr)).max()
        if m > worst:
            worst = m
    return worst
```

Add the imports at the top of the file (after the numpy import):

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
```

Replace the `print('Correctness check passed!')` block in `main()` with the defense check plus plots:

```python
    print('Correctness check passed!')

    # --- Defense: CPA-style correlation of traces with HW(s) ---
    challenges_flat = np.array([ch.ravel() for ch in challenges])
    hws = np.array([[hw(v) for v in row] for row in challenges_flat])

    traces_unmasked = sim_unmasked.get_traces()
    traces_masked = sim_masked.get_traces()

    corr_unmasked = max_corr_vs_secret(traces_unmasked, hws)
    corr_masked = max_corr_vs_secret(traces_masked, hws)
    print('Unmasked max |corr| vs HW(s): {:.3f}'.format(corr_unmasked))
    print('Masked   max |corr| vs HW(s): {:.3f}'.format(corr_masked))

    assert corr_unmasked > 0.4, 'Unmasked version should leak the secret!'
    assert corr_masked < corr_unmasked * 0.5, 'Masked version still leaks!'
    assert corr_masked < 0.5, 'Masked max correlation too high!'

    # --- Plots: the two traces must differ ---
    x = np.arange(traces_unmasked.shape[1])
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    axes[0].plot(x, traces_unmasked[0], color='C0')
    axes[0].set_title('Unmasked KyberNTT - power trace (trace 0)')
    axes[1].plot(x, traces_masked[0], color='C1')
    axes[1].set_title('Masked KyberNTT - power trace (trace 0)')
    for ax in axes:
        ax.set_xlabel('instruction index')
        ax.set_ylabel('power')
    fig.tight_layout()
    fig.savefig('traces_comparison.png', dpi=120)

    print('Saved plots to traces_comparison.png')
    print('All checks passed!')
```

- [ ] **Step 3: Run the full verification**

```bash
.venv/bin/python3 verify_masking.py
```

Expected: `Correctness check passed!`, then `Unmasked max |corr| vs HW(s): <high, > 0.4>` and `Masked max |corr| vs HW(s): <low, < 0.5>` with `Masked < 0.5 * Unmasked`, then `Saved plots to traces_comparison.png` and `All checks passed!`.

- [ ] **Step 4: Verify the plot file was written**

```bash
ls -la traces_comparison.png
```

Expected: a non-empty PNG file exists.

- [ ] **Step 5: Commit**

```bash
git add verify_masking.py traces_comparison.png
git commit -m "Add defense (leakage) verification and trace comparison plot"
```

---

## Notes

- If the Task 4 correctness assertion fails, the bug is almost certainly in the mod-q share split or the recombination in `project.c`. Re-check that `s0` and `s1` stay in `[0, q]` before the NTTs and that recombination uses `barrett_reduce` after adding the two (already reduced) NTT outputs.
- If the Task 5 masked correlation is too high, confirm the split happens entirely before `starttrigger()` and that printing happens after `endtrigger()`.
- The `verify_masking.py` script must be run from the repo root (it relies on `repository='elmo/projects'` and the source-tree `elmo` package).
