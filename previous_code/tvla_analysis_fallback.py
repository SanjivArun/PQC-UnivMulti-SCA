#!/usr/bin/env python3
"""
tvla_analysis.py -- TVLA (Test Vector Leakage Assessment) study of the KyberNTT
masked vs unmasked ELMO simulations, in the "fixed vs random" design.

For each version of Kyber (unmasked ``KyberNTTSimulation`` and masked
``KyberNTTMaskedSimulation``) this script:

1. Runs the simulation with ``--nb-fixed`` traces built from a *constant*
   (fixed) secret followed by ``--nb-random`` traces built from *random*
   secrets (must be equal, a balanced fixed-vs-random TVLA).
2. Loads the produced power traces.
3. Computes the **ground-truth** leakage count using ELMO's *built-in*
   fixed-vs-random t-test (run with ``-fvr <N>``), which reports the number of
   instruction points with ``|t| > 4.5``.
4. Computes the leakage count with our own **univariate Welch's t-test**
   (``tvla.count_univariate_leakages``).
5. Computes the leakage count with our own **multivariate Hotelling T^2 test**
    over non-overlapping windows of 50 points (``tvla.count_multivariate_leakages``).
6. Prints a side-by-side comparison table.

Expected outcome (sound 2-share masking):
    * unmasked:  large ground-truth / univariate counts, some leaking windows.
    * masked:    ~0 everywhere (first-order masking defeats fixed-vs-random TVLA).

Run from the repo root (the simulations share one ELMO output directory, so each
version's traces are read immediately after its own ``run()``):

    .venv/bin/python3 tvla_analysis.py            # 256 fixed + 256 random / version
    .venv/bin/python3 tvla_analysis.py --nb-fixed 128 --nb-random 128   # quick run
    .venv/bin/python3 tvla_analysis.py --verbose  # print ELMO progress/output

Full options: ``python3 tvla_analysis.py --help``.
"""

import argparse
import os
import re
import subprocess

import numpy as np

from elmo import get_simulation

import tvla

# --------------------------------------------------------------------------
# Paths / constants
# --------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
ELMO_TOOL_DIR = os.path.join(REPO_ROOT, "elmo", "elmo-tool")

UNMASKED_BIN = os.path.abspath(
    os.path.join(REPO_ROOT, "elmo", "projects", "Examples", "KyberNTT", "project.bin"))
MASKED_BIN = os.path.abspath(
    os.path.join(REPO_ROOT, "elmo", "projects", "Examples", "KyberNTTMasked", "project.bin"))

#: Per-version definition: (label, simulation class name, binary path).
VERSIONS = [
    ("unmasked", "KyberNTTSimulation", UNMASKED_BIN),
    ("masked", "KyberNTTMaskedSimulation", MASKED_BIN),
]


# --------------------------------------------------------------------------
# Challenge generation
# --------------------------------------------------------------------------

def make_challenges(nb_fixed, nb_random, fixed_value=0, k=2, n=256):
    """Build the fixed-vs-random challenge list.

    Returns ``[fixed] * nb_fixed + [random] * nb_random`` where every fixed
    challenge is the same constant matrix and every random challenge is drawn
    from Kyber's coefficient distribution ``{-2, -1, 0, 1, 2}``. Seeded for
    reproducibility (same fixed seed as the other repo scripts).
    """
    rng = np.random.RandomState(0)
    fixed = np.full((k, n), fixed_value, dtype=np.int64)
    challenges = [fixed.copy() for _ in range(nb_fixed)]
    for _ in range(nb_random):
        challenges.append(rng.choice(
            [-2, -1, 0, 1, 2], (k, n), p=[1 / 16, 4 / 16, 6 / 16, 4 / 16, 1 / 16]))
    return challenges


# --------------------------------------------------------------------------
# ELMO simulation + trace loading
# --------------------------------------------------------------------------

def clean_output_traces():
    """Delete stale per-trace outputs from the shared ELMO output directory.

    Both simulations write to the *same* ``elmo/elmo-tool/output/`` directory,
    so traces from the previous version (or from an earlier run) must be removed
    before the next simulation: a stale trace that is not overwritten by ELMO
    would otherwise be silently mixed into the new trace set.
    """
    for folder in ("traces", "nonprofiledindexes", "asmoutput"):
        path = os.path.join(ELMO_TOOL_DIR, "output", folder)
        if os.path.isdir(path):
            for name in os.listdir(path):
                os.remove(os.path.join(path, name))


def run_simulation(classname, challenges, verbose=False):
    """Run one simulation with the given challenges; return the simulation object."""
    clean_output_traces()
    Simulation = get_simulation(classname, repository="elmo/projects")
    simulation = Simulation()
    simulation.set_challenges(challenges)
    res = simulation.run()
    if res.get("error"):
        raise RuntimeError("ELMO run failed ({}): {}".format(classname, res["error"]))
    if verbose:
        print("  ran {}: {} traces x {} instructions".format(
            classname, res["nb_traces"], res["nb_instructions"]))
    return simulation


def load_traces(filenames):
    """Load trace files (one float per line) into a float32 ``(n, p)`` array.

    Loads directly from the trace files (bypassing ``get_traces``' in-memory
    list of Python floats) so the large masked set stays ~345 MB.
    """
    lengths = {}
    for fn in filenames:
        with open(fn) as fh:
            lengths.setdefault(fh.read().count("\n"), fn)
    if len(lengths) != 1:
        raise RuntimeError(
            "inconsistent trace lengths ({}) -- the trace directory was not "
            "cleanly rewritten by the simulation; re-run with a clean output "
            "directory".format(sorted(lengths.items())[:5]))
    length, first_fn = next(iter(lengths.items()))
    out = np.empty((len(filenames), length), dtype=np.float32)
    for i, fn in enumerate(filenames):
        out[i] = np.array(open(fn).read().split(), dtype=np.float32)
    return out


def run_elmo_fvr(binary_path, nb_fixed, verbose=False):
    """Ground truth: run ELMO's *built-in* fixed-vs-random t-test on the traces
    that are already in ``elmo/elmo-tool/output/traces/`` (traces 1..N are taken
    as "fixed", traces N+1..2N as "random").

    ``elmo <bin> -fvr N`` only *analyses* the existing trace files (it returns
    before emulating), so it must be run right after the simulation that
    produced them. Returns the parsed count of leaking instructions/cycles
    (points with |t| > 4.5).
    """
    cmd = '{} "{}" -fvr {}'.format(os.path.join(ELMO_TOOL_DIR, "elmo"), binary_path, nb_fixed)
    proc = subprocess.run(cmd, shell=True, cwd=ELMO_TOOL_DIR,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("ELMO -fvr failed (rc={}): {}".format(
            proc.returncode, proc.stderr[-500:]))
    if verbose:
        for line in proc.stdout.splitlines():
            if "LEAKY" in line or "first order" in line or "instructions/cylces" in line:
                print("  elmo: " + line)
    match = re.search(
        r"first order fixed vs random fail instructions/cycles\s+(\d+)", proc.stdout)
    if match is None:
        raise RuntimeError(
            "Could not parse ELMO -fvr output. stdout tail: {}".format(proc.stdout[-500:]))
    return int(match.group(1))


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def print_table(results, args):
    print()
    print("=" * 78)
    print("TVLA leakage detection -- fixed (first N traces) vs random (next N)")
    print("=" * 78)
    print("  traces per group : {}".format(args.nb_fixed))
    print("  univariate       : Welch |t| > {:.2f}".format(tvla.DEFAULT_TVLA_THRESHOLD))
    print("  multivariate     : Hotelling T^2, window={} pts, alpha={:.3e}".format(
           args.window_size, tvla.DEFAULT_ALPHA))
    print("-" * 78)
    header = "{:<12} {:>14} {:>13} {:>15}".format(
        "version", "ground truth", "univariate", "multivariate")
    subheader = "{:<12} {:>14} {:>13} {:>15}".format(
        "", "(ELMO -fvr)", "(Welch t)", "(Hotelling T2)")
    print(header)
    print(subheader)
    print("-" * 78)
    for name, r in results.items():
        print("{:<12} {:>14} {:>13} {:>15}".format(
            name, r["ground_truth"], r["univariate"], r["multivariate"]))
    print("=" * 78)
    for name, r in results.items():
        print("  {}: {} singular-covariance windows (of {} windows) "
              "handled with the pseudo-inverse".format(
                  name, r["n_singular"], r["n_windows"]))
        if r["n_singular"]:
            print("    (many singular windows are expected for the unmasked version, "
                  "whose 'fixed' traces are all identical)")


def _make_summary_figure(results, args):
    """Save a summary figure that is readable without zooming or technical background.

    Panel A (top): bar chart showing what fraction of trace points were flagged.
    Panel B (bottom): bar chart showing what fraction of windows were flagged.
    Each bar goes from 0 at the bottom; 0-height bars appear as a line on the axis.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=False)

    versions = ["unmasked", "masked"]
    colors = ["#c0392b", "#27ae60"]

    # ---- Panel A: leaked points ----
    ax = axes[0]
    for i, name in enumerate(versions):
        r = results[name]
        n_total = r["t_stats"].shape[0]
        leaked = r["univariate"]
        pct = leaked / n_total * 100 if n_total > 0 else 0
        yval = leaked
        ax.bar(name, yval, color=colors[i], width=0.5, edgecolor="white", linewidth=0.5)
        ax.text(i, yval / 2,
                "{:,} of {:,}\n({:.1f}%)\n-- LEAKING --".format(
                    int(leaked), n_total, pct),
                ha="center", va="center", fontsize=8, color="white", fontweight="bold")
    ax.set_ylim(0, None)
    ax.set_ylabel("Leaked trace points")
    ax.set_title("How many individual trace points leaked the secret? (univariate test)")

    # ---- Panel B: leaked windows ----
    ax = axes[1]
    max_windows = max(r["n_windows"] for r in results.values())
    for i, name in enumerate(versions):
        r = results[name]
        leaked = r["multivariate"]
        total_w = r["n_windows"]
        pct = leaked / total_w * 100 if total_w > 0 else 0
        ax.bar(name, leaked, color=colors[i], width=0.5, edgecolor="white", linewidth=0.5)
        ax.text(i, leaked / 2 if leaked > 0 else max_windows * 0.05,
                "{:,} of {:,}\n({:.1f}%)\n-- {} --".format(
                    int(leaked), total_w, pct, "LEAKING" if leaked > 0 else "SAFE"),
                ha="center", va="center", fontsize=8,
                color="white" if leaked > 0 else "black", fontweight="bold")
    ax.set_ylim(0, None)
    ax.set_ylabel("Leaked windows")
    ax.set_xlabel("")
    ax.set_title("How many 50-point windows leaked? (per-window multivariate test)")

    fig.suptitle(
        "KyberNTT TVLA -- side-by-side comparison: can it steal your Kyber key?",
        fontsize=11, fontweight="bold", y=0.96)

    ax.set_xticks(np.arange(len(versions)))

    fig.tight_layout(rect=[0, 0, 0.97, 0.93])
    return fig


def make_plots(results, args):
    """Save the technical detail plot (|t| per point + T² per window)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 7), sharex=False)
    threshold = tvla.DEFAULT_TVLA_THRESHOLD

    for col, (name, r) in enumerate(results.items()):
        # --- |t| per point ---
        ax_t = axes[0, col]
        ax_t.plot(np.arange(r["t_stats"].shape[0]), np.abs(r["t_stats"]),
                  color="C0", lw=0.5)
        ax_t.axhline(threshold, color="r", ls="--", lw=1,
                     label="Threshold = {}".format(threshold))
        ax_t.set_title("{}: {} |t|>{} points".format(name, r["univariate"], threshold))
        ax_t.set_ylabel("|t|")
        ax_t.set_yscale("symlog")
        ax_t.legend(fontsize=7)

        # --- T² per window ---
        ax_t2 = axes[1, col]
        # colour-coded: green below threshold, red above
        t2_vals = r["t2"]
        f_crit = tvla.f_critical(50, args.nb_fixed + args.nb_fixed - 2,
                                 tvla.DEFAULT_ALPHA)
        below = t2_vals < f_crit
        above = ~below
        max_t2 = max(t2_vals.max(), f_crit) if t2_vals.size > 0 else f_crit

        ax_t2.scatter(np.where(below)[0], t2_vals[below],
                      color="#2ecc71", s=1.5, edgecolors="none", label="Below threshold")
        ax_t2.scatter(np.where(above)[0], t2_vals[above],
                      color="#e74c3c", s=1.5, edgecolors="none", label="Above threshold")
        ax_t2.axhline(f_crit, color="gray", ls="--", lw=0.8,
                      label="T² threshold = {:.1f}".format(f_crit))
        ax_t2.set_title("{}: {} leaking windows".format(name, r["multivariate"]))
        ax_t2.set_ylabel("T²")
        ax_t2.set_xlabel("window index")
        ax_t2.set_ylim(0, max_t2 * 1.1)
        ax_t2.legend(fontsize=7)

    fig.tight_layout()
    out = os.path.join(REPO_ROOT, "tvla_results.png")
    fig.savefig(out, dpi=120)
    print("Saved plots to {}".format(out))

    # Save a simplified, annotated summary figure
    summary = _make_summary_figure(results, args)
    summary_out = os.path.join(REPO_ROOT, "tvla_summary.png")
    summary.savefig(summary_out, dpi=120)
    plt.close(summary)
    print("Saved summary to {}".format(summary_out))


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fixed-vs-random TVLA leakage study of masked vs unmasked KyberNTT.")
    parser.add_argument("--nb-fixed", type=int, default=256,
                        help="number of traces with the fixed (constant) secret (default 256)")
    parser.add_argument("--nb-random", type=int, default=256,
                        help="number of traces with random secrets (default 256)")
    parser.add_argument("--window-size", type=int, default=50,
                        help="multivariate window size in points (default 50)")
    parser.add_argument("--fixed-value", type=int, default=0,
                        help="constant secret value for the fixed traces (default 0)")
    parser.add_argument("--verbose", action="store_true",
                        help="print ELMO progress and -fvr output")
    parser.add_argument("--save-stats", action="store_true",
                        help="save per-version t-stats / T^2 to tvla_stats.npz")
    parser.add_argument("--no-plots", action="store_true",
                        help="do not save tvla_results.png")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.nb_fixed != args.nb_random:
        raise SystemExit(
            "--nb-fixed must equal --nb-random (balanced fixed-vs-random TVLA).")

    challenges = make_challenges(args.nb_fixed, args.nb_random, args.fixed_value)

    results = {}
    for label, classname, bin_path in VERSIONS:
        print("=== {} ({}) ===".format(label, classname))
        sim = run_simulation(classname, challenges, args.verbose)
        filenames = sim.get_results_filenames()
        print("  loading {} traces...".format(len(filenames)))
        traces = load_traces(filenames)

        print("  running ELMO -fvr (ground truth)...")
        ground_truth = run_elmo_fvr(bin_path, args.nb_fixed, args.verbose)

        group0 = traces[: args.nb_fixed]
        group1 = traces[args.nb_fixed: args.nb_fixed + args.nb_fixed]
        univariate, t_stats = tvla.count_univariate_leakages(group0, group1)
        multivariate, t2, n_singular, _f_crit = tvla.count_multivariate_leakages(
            group0, group1, args.window_size)

        results[label] = {
            "ground_truth": ground_truth,
            "univariate": univariate,
            "multivariate": multivariate,
            "t_stats": t_stats,
            "t2": t2,
            "n_singular": n_singular,
            "n_windows": t2.shape[0],
        }
        print("  ground truth: {} leaky points | univariate: {} | multivariate: {}"
              .format(ground_truth, univariate, multivariate))
        del traces  # free ~345 MB before running the next version

    print_table(results, args)

    if args.save_stats:
        np.savez(os.path.join(REPO_ROOT, "tvla_stats.npz"),
                 unmasked_t=results["unmasked"]["t_stats"],
                 masked_t=results["masked"]["t_stats"],
                 unmasked_t2=results["unmasked"]["t2"],
                 masked_t2=results["masked"]["t2"])
    if not args.no_plots:
        make_plots(results, args)


if __name__ == "__main__":
    main()
