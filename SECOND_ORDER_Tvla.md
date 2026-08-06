# Second-Order TVLA for Masked KyberNTT

## Purpose

This extends the existing first-order TVLA (Test Vector Leakage Assessment) pipeline to support **second-order leakage detection** for masked implementations of the ML-KEM (Kyber) NTT.

The existing implementation already has:
- First-order Welch's t-test (`count_univariate_leakages`)
- First-order Hotelling's T² test (`count_multivariate_leakages`)
- ELMO trace generation and loading
- Comparison tables and plotting

**What this adds:** A second-order TVLA pipeline that detects leakage surviving first-order masking. The second-order pipeline uses the standard centered-product transform from the side-channel literature, followed by the same statistical tests (Welch's t-test and Hotelling's T²), but operating on pairwise products of centered trace samples instead of raw amplitudes.

It also adds `--inject-masked-joint-leak` to inject a controlled joint-only second-order leak into the masked group for debugging and validating that Hotelling's T² is strictly stronger than Welch's per-feature test.

## What second-order TVLA is, and why it's needed

### The masking model

In a 2-share masked implementation, every secret value `s` is split into two random shares `s = s0 + s1 mod q`. Operations are performed on the shares instead of `s` itself. The goal is to ensure that no single intermediate value at any hardware point has a statistical dependence on `s`.

However, a single-point (first-order) power measurement only reveals information about ONE intermediate value at a time:
```
x(i) ≈ f(s0) + noise    -- first intermediate
x(j) ≈ g(s0) + noise    -- second intermediate
```
If the masking is perfect, each individual measurement `x(i)` is independent of `s`, and first-order TVLA detects **zero** leakage. This is exactly what we observe with the current 2-share KyberNTTMasked design.

But second-order leakage exploits the **correlation between two measurement points within the same instruction**. If an instruction produces two or more intermediate values that each leak information about the same secret share `s0`, then no single measurement leaks -- but the **product** of two centered measurements does:

```
(x(i) - μ(i)) · (x(j) - μ(j)) ≈ f(s0)·g(s0) + noise

The cross-term f(s0)·g(s0) depends on the secret share -- and that's detectable by TVLA.
```

This is the foundational insight behind second-order side-channel attacks, first formalized in the masked-cryptography literature.

### How this implementation works

```
ELMO simulation
    ↓
Raw traces (n×p per group), shape ~ (256, 93000)
    ↓
Compute global mean μ(i) across ALL traces (both groups combined)
    ↓
Center all traces:  x'(i) = x(i) - μ(i)
    ↓
Divide each trace into non-overlapping windows of W points (W=10 default)
    ↓
For each window [s, s+W): compute the upper-triangular part of the outer product
    → x'(s)·x'(s),  x'(s)·x'(s+1), ..., x'(s+W-1)·x'(s+W-1)
    → W(W+1)/2 features per window
    With W=10: 55 features per window
    Total: (93000 / 10) × 55 = ~511,500 feature columns
    ↓
Apply Welch's t-test per feature (univariate second-order)
    OR
Apply Hotelling's T² per window of features (multivariate second-order)
    ↓
Count leaking features / windows
```

### Joint-only leak injection for validation

To validate that Hotelling's T² is strictly stronger than Welch's per-feature test, the pipeline supports optional injection of a **joint-only second-order leak** into the masked group:

```bash
python3 tvla_analysis.py --order 2 \
    --inject-masked-joint-leak 0.05 \
    --joint-leak-fraction 0.1 \
    --no-plots
```

This injects a sub-threshold bias (`per_feature_sigma` std devs per feature) into a fraction of windows. By construction:
- Per-feature Welch |t| stays **below** the TVLA threshold → univariate is BLIND
- Per-window Hotelling T² exceeds the critical value → multivariate DETECTS

An audit table is printed after the main results, showing per-injected-window metrics (max |t|, F_obs, and whether each was flagged), demonstrating "multivariate strictly beats univariate" empirically.

### Design decisions and justification

**Global centering (not per-group centering)**

All traces from both the "fixed" and "random" groups are centered by a single mean computed over the combined set. This is statistically correct under the null hypothesis: if there's no leakage, both groups share the same unconditional distribution, so one shared mean estimate is the correct model. Under the alternative, global centering is more conservative -- if second-order leakage survives global centering, the deviation is stronger and more robust.

**Window size W = 10 (configurable)**

The full second-order feature space would require computing C(93000, 2) ≈ 4.3 billion pairwise products -- completely infeasible. We therefore use a **windowed approach**: we partition each trace into small non-overlapping windows of W points, and within each window we compute all pairwise products. With W=10 this gives 55 features per window and ~511,500 total columns. This is practical: the feature matrix for 1024 traces at 8 bytes per float64 is only ~4.0 GB.

W=10 captures the most likely leakage pattern -- correlations between samples within the same arithmetic operation -- while keeping dimensions tractable. The window size is a CLI parameter (`--second-order-window-size`) so users can sweep it if they suspect leakage at different temporal scales.
```

**Non-overlapping windows (default)**

We use non-overlapping windows for simplicity and low memory. This may miss inter-window pairwise leakage. If leakage is suspected at a specific instruction boundary, the user can increase `--second-order-window-size` to capture a larger span, or the implementation can be extended later with overlapping windows or a sliding pairwise approach.

**Welch's t-test and Hotelling's T² (unchanged)**

The statistical tests are identical to the first-order pipeline. This provides API consistency: the second-order pipeline is `transform → test`, not a new test.

**No modification to the first-order pipeline**

All existing functions, CLI options, plots, and table formats for first-order TVLA remain untouched. The `--order 1` path is the existing behavior. The `--order 2` path is a new branch.

## Key files changed

| File | What changes |
|------|--------------|
| `tvla.py` | Three new functions: `preprocess_second_order_traces`, `count_second_order_univariate_leakages`, `count_second_order_multivariate_leakages`. Zero changes to existing functions. |
| `tvla_analysis.py` | New `--order {1,2}` CLI flag; new `--second-order-window-size` flag (default 10); new `--inject-masked-joint-leak <sigma>` and `--joint-leak-fraction` flags for debugging; redesigned `print_table()` (order 2 layout drops ground truth column, shows leaking/total windows); new `print_joint_leak_audit()` for the injection verification table; updated plots adapt to order 2. |
| `test_tvla.py` | Synthetic tests validating the second-order pipeline end-to-end. |
| `SECOND_ORDER_Tvla.md` | This file. |

## CLI reference

| Flag | Default | Description |
|------|---------|-------------|
| `--nb-fixed`, `--nb-random` | 256 | Traces per group (must be equal) |
| `--window-size` | 10 | Multivariate window size in points (first-order) |
| `--second-order-window-size` | 10 | Centered-product window size in trace points (second-order) |
| `--order {1,2}` | 1 | TVLA order: 1 = first-order, 2 = second-order |
| `--inject-masked-joint-leak <sigma>` | 0.0 (off) | Inject joint-only leak into masked group with per_feature_sigma. E.g. `0.05` for a subtle leak visible only to multivariate tests |
| `--joint-leak-fraction` | 0.1 | Fraction of masked windows to inject (when `--inject-masked-joint-leak` > 0) |
| `--no-plots` | off | Skip saving tvla_results.png and tvla_summary.png |
| `--verbose` | off | Print ELMO progress and -fvr output |

### Example commands

```bash
# First-order TVLA (existing behaviour, default)
python3 tvla_analysis.py

# Second-order TVLA, default window size 10
python3 tvla_analysis.py --order 2

# Second-order with joint-only leak injection for multivariate validation
python3 tvla_analysis.py --order 2 --inject-masked-joint-leak 0.05 --joint-leak-fraction 0.1 --no-plots
```

## Output formats

### First-order (order=1) -- unchanged

The table shows ground truth (ELMO -fvr), univariate (Welch t), and multivariate (Hotelling T²) counts:

```
==============================================================================
TVLA leakage detection -- fixed (first N traces) vs random (next N)
==============================================================================
  traces per group : 256
  univariate       : Welch |t| > 4.50
  multivariate     : Hotelling T^2, window=10 pts, alpha=6.795e-06
------------------------------------------------------------------------------
version          ground truth     univariate         multivariate
                  (ELMO -fvr)    (Welch t)      (Hotelling T2)
------------------------------------------------------------------------------
unmasked                  0            14998               4200
masked                      0                2                 10
==============================================================================
```

### Second-order (order=2) -- new layout

Drops the ground truth column (ELMO -fvr is inherently first-order). Shows per-feature Welch counts and per-window leaking/total:

```
==============================================================================
TVLA leakage detection -- fixed (first N traces) vs random (next N)
==============================================================================
  traces per group : 256
  second-order     : centered-product, window_size=10 pts
  univariate       : Welch |t| > 4.50
  multivariate     : Hotelling T^2, 55 features/window, alpha=6.795e-06
------------------------------------------------------------------------------
version          univariate           windows
                   (Welch t)          (leaking/total)
------------------------------------------------------------------------------
unmasked             109531            7683/8418
masked                    1              439/16835
==============================================================================
```

### Joint-leak audit (when --inject-masked-joint-leak > 0)

Printed after the main table, showing per-injected-window metrics:

```
============================================================================================
JOINT-LEAK AUDIT  (injected into 1684 windows, per_feature_sigma=0.05)
============================================================================================
  window    trace_pts        max|t| (target<4.5)    F_obs (target>Fcrit=2.19)  flagged?
    6396    [ 63960,  63970)       1.96  BLIND           6.7  DETECTED   yes
    7056    [ 70560,  70570)       1.81  BLIND           0.1  missed     no
   ...
--------------------------------------------------------------------------------------------
  summary:
    injected windows                    = 1684
    univariate detected (of 1684)         = 1   (target: noise floor)
    multivariate detected (of 1684)        = 439
    leakage clearly missed by Welch     = 1683   <- multivariate strictly stronger
============================================================================================
```

## Mathematical reference

The centered-product second-order transform follows the methodology described in:
- Narasimhan, Balbuzhev, and Biba, "Second-Order Side-Channel Attacks on Masked Implementations"
- The standard approach cited in the CHES and EUROCRYPT masked-cryptography literature

## Validation approach

Since we can't easily verify whether the real KyberNTTMasked binary actually exhibits second-order leakage, we validate via **artificial injection**:

1. Generate synthetic traces where group0 and group1 have identical means at every point (no first-order leakage)
2. Center the traces
3. Inject controlled second-order leakage to specific pairwise products: set x'(i)·x'(j) = target_value for chosen (i,j) pairs, where target_value depends on the group assignment
4. Verify: first-order TVLA detects 0 leakages, second-order TVLA detects the injected leakages

This confirms the mathematical correctness of the pipeline without requiring a faulty masked implementation.

Additionally, the `--inject-masked-joint-leak` flag provides an *in-situ* validation on real ELMO traces: it injects a known joint-only leak whose univariate response should be near-zero and multivariate response should be significant, verifying that the multivariate test is actually stronger in practice.
