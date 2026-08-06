# How to run the KyberNTT example

This document explains how to install `python-elmo` into the project venv and run the
KyberNTT simulation example (`elmo/projects/Examples/KyberNTT/`), the masked variant
(`elmo/projects/Examples/KyberNTTMasked/`), the masking verification script
(`verify_masking.py`), and the fixed-vs-random TVLA study (`tvla_analysis.py`). It also
explains what `setup.py` and `test.py` are for.

## Prerequisites

- Python 3 (use `python3`, not `python`)
- `gcc`, `make`, and `git` on `PATH` (the install step compiles the ELMO C tool)
- Internet access during installation (the ELMO tool is cloned from GitHub)

## 1. Activate the venv

```bash
source .venv/bin/activate
```

When the venv is active, `python3` and `python3 -m pip` already point at the venv.
If you do not activate it, prefix every command below with `.venv/bin/` instead.

## 2. Install the Python dependencies

```bash
python3 -m pip install numpy matplotlib
```

## 3. Install python-elmo (this also compiles the ELMO tool)

```bash
python3 -m pip install .
```

This is the key step. Besides installing the Python package into the venv, `setup.py`
runs a custom build step that:

1. clones the ELMO C tool (`https://github.com/sca-research/ELMO`) into `elmo/elmo-tool`,
2. patches `elmodefines.h` so the tool reads/writes `input.txt` (the file python-elmo uses),
3. runs `make` to compile the `elmo` simulator binary.

If you later change the C sources or switch machines, re-run this step.

## 4. Run the KyberNTT simulation

Run this script (or paste it into an interactive `python3` shell):

```python
from elmo import get_simulation

KyberNTTSimulation = get_simulation('KyberNTTSimulation')

simulation = KyberNTTSimulation()
challenges = simulation.get_random_challenges(10)
simulation.set_challenges(challenges)

simulation.run()            # launches the ELMO simulator
traces = simulation.get_traces()
print(traces.shape)         # (10, nb_instructions)
```

The class is defined in `elmo/projects/Examples/KyberNTT/projectclass.py`. It simulates
the NTT of Kyber512 (K=2, N=256) and, for each challenge (a 2x256 coefficient matrix),
records the power trace of the leaking C code in `project.c`.

Useful variants:

- Deterministic inputs: `simulation.get_test_challenges()` instead of random ones.
- `simulation.get_indexes_of(lambda instr: 'mul' in instr)` to get the indexes of the
  MUL instructions, and `simulation.get_traces(indexes)` to slice the traces.
- `simulation.run_online()` to run through the client-server path
  (`python3 -m elmo run-server` must be running).

## 5. (Optional) Quick end-to-end check

```bash
python3 test.py
```

This installs the ELMO tool if missing, then runs the full set of smoke tests,
including the KyberNTT simulation.

## 6. The masked variant: KyberNTTMasked

`elmo/projects/Examples/KyberNTTMasked/` is a side-channel-protection study of the same
NTT. It splits the secret `s` into two random shares `s0, s1` (with `s0+s1 = s`), NTTs
each share separately, and recombines the two NTT outputs mod q *after* the recorded
trace window closes, so the true NTT(s) is never present in the leaked power trace.

Source layout (same structure as the unmasked example, plus the masking):

- `project.c` — the "secretly leaking" main program. `starttrigger()`/`endtrigger()`
  delimit the recorded window. The share split happens before `starttrigger()`, both
  NTTs run inside the window, and the recombination (`barrett_reduce(s0+s1)`) runs
  *after* `endtrigger()`.
- `projectclass.py` — defines `KyberNTTMaskedSimulation`; reuses the challenge format
  of the unmasked example, and additionally records the random bytes ELMO drew for the
  share split (`randdata.txt`), so a run can be reproduced bit-for-bit.

Run it the same way as the unmasked example:

```python
from elmo import get_simulation

Sim = get_simulation('KyberNTTMaskedSimulation', repository='elmo/projects')
sim = Sim()
sim.set_challenges(sim.get_random_challenges(10))
sim.run()
print(sim.get_traces().shape)          # (10, nb_instructions) — note: longer than unmasked
print(sim.get_printed_data())
```

### Rebuilding the masked binary

If you change any C source, rebuild the binary (the prebuilt `project.bin` is committed
so the example runs without a cross-toolchain):

```bash
make -C elmo/projects/Examples/KyberNTTMasked clean
make -C elmo/projects/Examples/KyberNTTMasked
```

This uses the included prebuilt `arm-none-eabi` toolchain. On a machine without it, keep
using the committed `project.bin` and edit only the Python side.

## 7. Full masking verification: `verify_masking.py`

The end-to-end check that the masking is *both* correct and side-channel-safe. Run it
from the repo root:

```bash
.venv/bin/python3 verify_masking.py
```

What it does:

1. Runs both the unmasked and masked simulations with the same 256 random challenges
   (small secrets drawn from Kyber's distribution). Each simulation's printed data and
   traces are captured immediately after its own `run()`, because the two simulations
   share a single ELMO output directory that gets overwritten by each run.
2. **Correctness** — rebuilds the printed uint16 coefficients in big-endian order
   (`pd[2k] * 256 + pd[2k + 1]`, matching `print2bytes`) and asserts the masked NTT
   output equals the unmasked one modulo q=3329 for all 131072 coefficients.
3. **Defense** — a CPA-style check: Pearson correlation of every power-trace column
   against the Hamming weight of the secret, for both versions. It asserts the unmasked
   version leaks (max |corr| > 0.5) and the masked version does not
   (max |corr| < 0.5; in practice ~0.000, i.e. the noise floor).
   Columns with zero variance are skipped, since a constant column produces a spurious
   huge correlation.
4. **Plots** — saves `traces_comparison.png` (top: unmasked trace 0, bottom: masked
   trace 0) and prints `All checks passed!`.

Expected output (numbers may wobble slightly because ELMO reseeds its RNG from the clock
each run; the masked max |corr| tracks the null noise floor `~sqrt(2*ln(P)/N)`, so it
stays well below the 0.5 threshold):

```
Correctness: coefficients differ in 0 of 131072
Correctness check passed!
Unmasked max |corr| vs HW(s): 0.976
Masked   max |corr| vs HW(s): 0.000
Saved plots to traces_comparison.png
All checks passed!
```

### A note on the threshold

`max |corr|` for a *masked* trace is pure noise; for a given number of points `P` and
traces `N` it sits around `sqrt(2*ln(P)/N)`. That is why the script needs a few hundred
traces: with too few, the noise floor can exceed 0.5 and the (correct) masking would
still trip the defense assertion.

---

## 8. Fixed-vs-random TVLA: `tvla_analysis.py`

A fixed-vs-random TVLA (Test Vector Leakage Assessment) study comparing the masked and
unmasked KyberNTT. It supports **two orders of detection**:

### First-order (`--order 1`, the default)

For each version it runs the simulation with `--nb-fixed` traces of a
*constant* secret followed by `--nb-random` traces of *random* secrets, then reports:

| method | what it counts |
| --- | --- |
| **ground truth** | ELMO's *built-in* fixed-vs-random t-test (`elmo <bin> -fvr <N>`), i.e. points with `|t| > 4.5` |
| **univariate** | our own Welch t-test per point, `|t| > 4.5` (`tvla.py`) |
| **multivariate** | our own Hotelling T² per non-overlapping window of 50 points, tested at `alpha = 6.795e-06` (the two-sided `p` matching `|t| > 4.5`) |

The univariate count is checked against ELMO's ground truth (they match exactly, which
cross-validates `tvla.py`). Expected result: the unmasked version leaks thousands of points
(44,616 of 93,142 at 256 traces/group), the masked version leaks none (0 everywhere).

#### Run

```bash
.venv/bin/python3 tvla_analysis.py                               # first-order, 256 fixed + 256 random per version
.venv/bin/python3 tvla_analysis.py --nb-fixed 32 --nb-random 32  # quick run
.venv/bin/python3 tvla_analysis.py --verbose                     # print ELMO progress + -fvr output
```

#### Expected output (full run)

```
========================================================================
TVLA leakage detection -- fixed (first N traces) vs random (next N)
========================================================================
  traces per group : 256
  univariate       : Welch |t| > 4.50
  multivariate     : Hotelling T^2, window=50 pts, alpha=6.795e-06
------------------------------------------------------------------------------
version            ground truth    univariate          windows
                   (ELMO -fvr)     (Welch t)     (leaking/total)
------------------------------------------------------------------------------
unmasked                 44616         44616           206/1862
masked                       0             0             0/1683
========================================================================
```

### Second-order (`--order 2`)

After a first-order masking proof, the next question is: does the masking also defeat
_second-order_ attacks that look at products of nearby power samples (pairwise products
of centered samples from the same instruction)?

The second-order pipeline:

1. **Global centering**: compute one mean across ALL traces from both groups, subtract.
2. **Windowed centered-product transform**: for each window, compute all pairwise products
   of centered samples (upper triangular: `c[i]*c[j]` for `0 <= i <= j < W`). Each window
   of `W` points becomes `W*(W+1)/2` "second-order features".
3. **Welch t-test per feature** — univariate second-order detection.
4. **Hotelling T² per window** — multivariate second-order detection: treats each window
   of `W*(W+1)/2` features as a single multivariate block.

No ELMO `-fvr` ground truth is computed in second-order mode (it is inherently first-order
only). The masked version is expected to leak near-zero under both metrics (proving
protection against the centered-product attack too).

The unmasked version should produce the same leakage counts as the first-order pipeline in
univariate mode (the Welch t-test on centered products of the _same_ trace point `c[i]*c[i]`
is equivalent to the variance, which is exactly what the first-order Welch test compares
between groups), and potentially additional leakages in cross-terms `c[i]*c[j]` where `i ≠ j`
if nearby samples in the same window are jointly correlated with the secret.

#### Run

```bash
.venv/bin/python3 tvla_analysis.py --order 2                               # 2nd order, 256/256 traces, window size 10 (default)
.venv/bin/python3 tvla_analysis.py --order 2 --second-order-window-size 10 # explicit window size 10 (same as above)
.venv/bin/python3 tvla_analysis.py --order 2 --nb-fixed 64 --nb-random 64  # quick run
.venv/bin/python3 tvla_analysis.py --order 2 --second-order-window-size 5  # shorter windows (15 features/window)
```

#### Validating the detector with an injected leak

To prove the multivariate (and univariate) test actually fires on the real masked
traces, you can inject an *intentional* second-order leak into the masked random group
(`--order 2` only; the masked version normally reports `0/16,835`):

```bash
.venv/bin/python3 tvla_analysis.py --order 2 --inject-masked-leak 0.5 --leak-fraction 0.1
.venv/bin/python3 tvla_analysis.py --order 2 --inject-masked-leak 1.0 --leak-window-index 0  # leak exactly window 0
```

- `--inject-masked-leak <mag>` adds a bias of `<mag>` **standard deviations** to all 55
  centered-product features of the chosen window(s) in the masked *random* group (a
  one-sigma shift in `magnitude=1.0` is detectable regardless of the absolute power
  scale). That bias is a second-order leak, so the window column toggles from `0/N`
  to roughly `fraction×N / N`.
- `--leak-fraction <0..1>` leaks that fraction of windows (default `0.1`).
- `--leak-window-index <i>` leaks exactly one window (overrides `--leak-fraction`).
- A `NOTE:` line in the table marks the masked result as artificial.

Guards (fail fast, before the simulation runs): injection requires `--order 2`;
magnitude must be `> 0`; `--leak-fraction` must be in `(0, 1]`; a given
`--leak-window-index` must be within `[0, n_windows)` (checked once the trace length
is known).

#### Joint-only leak: where multivariate strictly beats univariate

The mean-shift above is detectable by *both* the univariate and multivariate tests,
so it validates the pipeline but does not separate them. To answer the research
question — *can a multivariate framework give stronger guarantees than univariate
TVLA?* — `--inject-masked-joint-leak` injects a **joint-only** leak: a per-feature bias
small enough that Welch stays blind (`|t| < 4.5` on every feature) but large enough that
the per-window Hotelling T² (Mahalanobis distance over all 55 features jointly) clears
its F critical value. This is the canonical "multivariate wins" construction.

```bash
.venv/bin/python3 tvla_analysis.py --order 2 --inject-masked-joint-leak 0.05 --joint-leak-fraction 0.1
.venv/bin/python3 tvla_analysis.py --order 2 --inject-masked-joint-leak 0.05 --joint-leak-window-index 42
.venv/bin/python3 tvla_analysis.py --order 2 --inject-masked-joint-leak 0.05 --joint-leak-fraction 0.1 --full-audit
```

- `--inject-masked-joint-leak <c>` injects a `c`-std per-feature bias into the masked
  random group. **Default `0.05`** — calibrated for `n = 256` traces/group so that
  per-feature Welch `|t| ≈ c·sqrt(n/2)` stays around 2.3 (safely below 4.5) while the
  multivariate Mahalanobis distance `D² ≈ features_per_window · c²` makes `F_obs` clear
  `F_crit` by a wide margin.
- `--joint-leak-fraction` (default `0.1`) / `--joint-leak-window-index` select the
  windows, same semantics as the mean-shift flags.
- `--full-audit` dumps every injected window in the audit (default: head + tail + summary).
- Mutually exclusive with `--inject-masked-leak` (mean-shift).

After the comparison table, a `JOINT-LEAK AUDIT` block reports, for each injected
window: its **location on the raw ELMO trace** in numpy index terms
`[w·W, (w+1)·W)`, the max per-feature Welch `|t|` (`BLIND` if < 4.5), and the per-window
Hotelling `T²` (`DETECTED` if ≥ `F_crit`). The summary row
`leakage clearly missed by Welch` is the empirical "yes" to the research question.

Phase 1 (synthetic, no ELMO) headline output:

```
JOINT-LEAK AUDIT  (injected into 10 windows, per_feature_sigma=0.05)
  window    trace_pts        max|t| (target<4.5)    T^2     (target>Fcrit=2.19)  flagged?
      26    [   260,    270)       3.83  BLIND         121.6  DETECTED   yes
      ...
  summary:
    injected windows                = 10
    univariate detected (of 10)     = 0   <- Welch BLIND to all
    multivariate detected (of 10)   = 10  <- Hotelling DETECTS all
    leakage clearly missed by Welch = 10  <- multivariate strictly stronger
```

Phase 2 (real masked KyberNTT) carries a caveat: every masked window has rank-deficient
covariance under the ELMO affine power model, so the pseudo-inverse collapses the
Mahalanobis distance and the empirical detection rate `M/K` may fall well below the
synthetic `10/10`. That honest number is itself a contribution — see `test_tvla.py`'s
joint-leak tests for the controlled statistical statement, and the Phase-2 run for the
real-trace calibration.

#### Expected output (full run)

```
========================================================================
TVLA leakage detection -- fixed (first N traces) vs random (next N)
========================================================================
  traces per group : 256
  second-order     : centered-product, window_size=10 pts
  univariate       : Welch |t| > 4.50
  multivariate     : Hotelling T^2, 55 features/window, alpha=6.795e-06
-----------------------------------------------------------------------------
version            univariate          windows
                    (Welch t)     (leaking/total)
------------------------------------------------------------------------------
unmasked              116742           8192/9314
masked                     1              0/16835
========================================================================
```

### CLI arguments

| argument | default | description |
| --- | --- | --- |
| `--nb-fixed` | 256 | number of fixed-secret traces per version |
| `--nb-random` | 256 | number of random-secret traces per version |
| `--window-size` | 50 | multivariate window size for first-order Hotelling T² |
| `--fixed-value` | 0 | constant secret value for fixed traces |
| `--order` | 1 | `1` = first-order (default), `2` = second-order centered-product |
| `--second-order-window-size` | 10 | window size for second-order centered-product transform |
| `--inject-masked-leak` | off | inject an intentional second-order leak of this magnitude into the masked random group (order 2 validation; mean-shift, both tests see it) |
| `--leak-fraction` | 0.1 | fraction of windows to leak with `--inject-masked-leak` |
| `--leak-window-index` | off | leak exactly this window index (overrides `--leak-fraction`) |
| `--inject-masked-joint-leak` | off | inject a JOINT-ONLY second-order leak of `c` std-devs/feature into the masked random group — Welch stays blind, Hotelling detects (the "multivariate strictly beats univariate" case) |
| `--joint-leak-fraction` | 0.1 | fraction of windows to leak with `--inject-masked-joint-leak` |
| `--joint-leak-window-index` | off | leak exactly this window index (overrides `--joint-leak-fraction`) |
| `--full-audit` | off | with `--inject-masked-joint-leak`, dump every injected window in the audit (default: head + tail + summary) |
| `--verbose` | off | print ELMO progress + `-fvr` output |
| `--save-stats` | off | dump raw t-stats / T² to `tvla_stats.npz` |
| `--no-plots` | off | skip saving `tvla_results.png` and `tvla_summary.png` |

### Files

- `tvla_analysis.py` — the CLI driver. Runs both simulations (cleaning the shared ELMO
  output directory first, because the two versions overwrite the same `output/traces/`),
  loads the traces, runs analysis at the selected order, and prints the comparison table.
  Saves `tvla_results.png` (per-point `|t|` and per-window T² for both versions) and
  `tvla_summary.png` (annotated summary bars) unless `--no-plots`.
- `tvla.py` — pure-numpy statistics module, independent of ELMO/scipy:
  - **First-order**: `welch_t_statistic`, `count_univariate_leakages`,
    `hotelling_t2_per_window`, `count_multivariate_leakages`
  - **Second-order**: `preprocess_second_order_traces` (global center + windowed product),
    `count_second_order_univariate_leakages`, `count_second_order_multivariate_leakages`
  - **Leak injection (validation)**: `inject_second_order_leak` (mean-shift, both tests see it),
    `inject_second_order_joint_leak_with_indices` (joint-only: Welch blind, Hotelling detects;
    returns the injected window indices for the per-window audit)
  - **Helpers**: `_extract_upper_triangular_indices`, `f_critical`, `_betacf`,
    `_regularized_beta`, `_beta_quantile`
- `test_tvla.py` — 30 manual checks (Welch t-test edge cases, F/beta reference values,
  leakage counts on synthetic data, comprehensive second-order tests, mean-shift leak
  injection, and the joint-only leak demonstrating multivariate strictly beats univariate).
  Run with `.venv/bin/python3 test_tvla.py`.
  Run with `.venv/bin/python3 test_tvla.py`.

---

## What are setup.py and test.py for?

### setup.py

A standard setuptools installation script (package name `python-elmo` 0.1.0, depends on
`numpy`). Its special part is the `PostBuildCommand` class (subclass of setuptools'
`build_py`): after the Python package is built, its `run()` calls
`install_elmo_tool()`, which clones, patches, and `make`-compiles the ELMO C simulator
into `elmo/elmo-tool`. So a single `python3 -m pip install .` gives you both the Python
wrapper and the compiled simulator that the wrapper shells out to at runtime.

### test.py

A manual smoke test (not a pytest/CI suite) meant to be run from the repo root to
verify an installation works end-to-end. It installs the ELMO tool if it is missing
(reusing `install_elmo_tool` from `setup.py`), then runs four checks and prints green
`Success!` messages:

1. **Manage a fresh simulation** — creates a new simulation project via
   `create_simulation`, finds it with `search_simulations`/`get_simulation`, then cleans up.
2. **Use the ELMO Engine** — computes power for 256 instruction points and asserts the
   result shape is `(256,)`.
3. **Use a real simulation** — runs the KyberNTT example with 10 random challenges,
   asserts the run has no error and produced 10 traces, then exercises the analysis
   helpers (`get_indexes_of`, `get_asmtrace`, `get_traces`, `get_printed_data`).
4. **Use ELMO by running online** — repeats test 3 through the client-server path
   (`launch_executor` / `run_online`).

If all four pass it prints `All seems fine!`. It exists so you can quickly confirm the
whole stack (Python wrapper + compiled ELMO tool + server mode) works after installing.
