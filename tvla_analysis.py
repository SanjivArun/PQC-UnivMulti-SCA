#!/usr/bin/env python3
"""
tvla_analysis.py -- TVLA (Test Vector Leakage Assessment) study of the Kyber
masked vs unmasked ELMO simulations, in the "fixed vs random" design.

For each version of Kyber (unmasked ``KyberKEMSimulation`` and masked
``KyberKEMMaskedSimulation``) this script:

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
    over non-overlapping windows of ``--window-size`` points (``tvla.count_multivariate_leakages``).
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
import trace_plot

# --------------------------------------------------------------------------
# Paths / constants
# --------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
ELMO_TOOL_DIR = os.path.join(REPO_ROOT, "elmo", "elmo-tool")

UNMASKED_BIN = os.path.abspath(
    os.path.join(REPO_ROOT, "elmo", "projects", "Examples", "Kyber", "project.bin"))
MASKED_BIN = os.path.abspath(
    os.path.join(REPO_ROOT, "elmo", "projects", "Examples", "KyberMasked", "project.bin"))

#: Per-version definition: (label, simulation class name, binary path).
VERSIONS = [
    ("unmasked", "KyberKEMSimulation", UNMASKED_BIN),
    ("masked", "KyberKEMMaskedSimulation", MASKED_BIN),
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
    err = res.get("error") or ""
    # ELMO emits benign ARM-emu stack diagnostics ("push {lr} ... popped 0x")
    # on stderr; they don't stop trace generation, so ignore unless traces
    # failed to be produced.
    real_err = "\n".join(
        line for line in err.splitlines() if "push {lr}" not in line
    ).strip()
    if real_err or not res.get("nb_traces"):
        raise RuntimeError("ELMO run failed ({}): {}".format(
            classname, err or "(no traces produced)"))
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
    if args.order == 1:
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
    else:
        W = args.second_order_window_size
        fpw = W * (W + 1) // 2
        print("=" * 78)
        print("TVLA leakage detection -- fixed (first N traces) vs random (next N)")
        print("=" * 78)
        print("  traces per group : {}".format(args.nb_fixed))
        print("  second-order     : centered-product, window_size={} pts".format(W))
        print("  univariate       : Welch |t| > {:.2f}".format(tvla.DEFAULT_TVLA_THRESHOLD))
        print("  multivariate     : Hotelling T^2, {} features/window, alpha={:.3e}".format(
               fpw, tvla.DEFAULT_ALPHA))
        print("-" * 78)
        header = "{:<12} {:>14} {:>18}".format("version", "univariate", "windows")
        subheader = "{:<12} {:>14} {:>18}".format("", "(Welch t)", "(leaking/total)")
        print(header)
        print(subheader)
        print("-" * 78)
        for name, r in results.items():
            print("{:<12} {:>14} {:>18}".format(
                name,
                "{:,}".format(r["univariate"]),
                "{:,}/{:,}".format(r["multivariate"], r["n_windows"])))
        print("=" * 78)

    print("-" * 78)
    for name, r in results.items():
        print("  {}: {} singular-covariance windows (of {} windows) "
              "handled with the pseudo-inverse".format(
                  name, r["n_singular"], r["n_windows"]))
        if r["n_singular"]:
            print("    (many singular windows are expected for the unmasked version, "
                  "whose 'fixed' traces are all identical)")


def print_joint_leak_audit(masked_result, inject_meta, args):
    """Print the joint-leak audit table."""
    W = args.second_order_window_size
    features_per_window = inject_meta["features_per_window"]
    n_injected = inject_meta["n_injected"]
    sigma = inject_meta["sigma"]
    inject_indices = inject_meta["inject_indices"]
    f_crit = inject_meta.get("f_crit")

    print()
    print("=" * 120)
    print("JOINT-LEAK AUDIT  (injected into {} windows, per_feature_sigma={:.2f})".format(n_injected, sigma))
    print("=" * 120)
    print("  {:>6}  {:>16}  {:>17}  {:>14}  {:>8}".format(
        "window", "trace_pts", "max|t| (target<4.5)", "F_obs (target>Fcrit={:.1f})".format(f_crit), "flagged?"))
    print("  " + "-" * 106)

    # Compute per-window metrics
    window_metrics = []
    t_stats = masked_result["t_stats"]
    t2 = masked_result["t2"]

    for w in inject_indices:
        w = int(w)
        w_start = w * features_per_window
        w_end = w_start + features_per_window
        t_start = w * W
        t_end = t_start + W
        max_abs_t = float(np.abs(t_stats[w_start:w_end]).max())
        nu = 2 * args.nb_fixed - features_per_window - 1
        scale = (2 * args.nb_fixed - 2) * features_per_window / nu if nu > 0 else 1.0
        f_obs = float(t2[w] / scale) if t2.shape[0] > w else 0.0
        flagged = t2[w] > f_crit if t2.shape[0] > w else False
        window_metrics.append({
            "window": w,
            "trace_range": "[{}, {})".format(t_start, t_end),
            "max_abs_t": max_abs_t,
            "f_obs": f_obs,
            "flagged": flagged,
        })

    # Sort: detected first (descending F_obs), then missed (descending |t|)
    detected = sorted(
        [m for m in window_metrics if m["flagged"]],
        key=lambda m: -m["f_obs"])
    missed = sorted(
        [m for m in window_metrics if not m["flagged"]],
        key=lambda m: -m["max_abs_t"])

    for m in (detected[:5] + missed[:5]):
        status = "DETECTED" if m["flagged"] else "missed"
        blind = "BLIND" if m["max_abs_t"] < 4.5 else ""
        print("  {:>6}  {:>16}  {:>6.2f}  {:>6} {:>8.1f}  {:>10}  {}".format(
            m["window"], m["trace_range"],
            m["max_abs_t"], blind,
            m["f_obs"], status,
            "yes" if m["flagged"] else "no"))

    if len(detected) + len(missed) > 10:
        print("  ...")

    multi_detected = sum(1 for m in window_metrics if m["flagged"])

    # Count leaking features ONLY within injected windows (for the univariate test)
    inject_threshold = tvla.DEFAULT_TVLA_THRESHOLD
    t_stats = masked_result["t_stats"]
    n_injected_features = n_injected * features_per_window
    uni_detected_injected = 0
    for w in inject_indices:
        w = int(w)
        w_start = w * features_per_window
        w_end = w_start + features_per_window
        leaked_in_window = np.count_nonzero(np.abs(t_stats[w_start:w_end]) > inject_threshold)
        uni_detected_injected += leaked_in_window
    print("-" * 120)
    print("  summary:")
    print("    injected windows                  = {}".format(n_injected))
    print("    injected features                 = {}".format(n_injected_features))
    print("    univariate detected (of {})       = {}   (target: noise floor)".format(n_injected_features, uni_detected_injected))
    print("    multivariate detected (of {})     = {} windows".format(n_injected, multi_detected))
    print("=" * 120)


def _make_summary_figure(results, args):
    """Save a summary figure that is readable without zooming or technical background.

    Panel A (top): bar chart of univariate (simple test) leaks per version.
    Panel B (bottom): bar chart of multivariate (advanced test) leaks per version.
    Each bar goes from 0 at the bottom; 0-height bars appear as a line on the axis.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=False)

    versions = ["unmasked", "masked"]
    colors = ["#c0392b", "#27ae60"]

    order_label = "second-order" if args.order == 2 else "first-order"

    # ---- Panel A: leaked features/points ----
    ax = axes[0]
    for i, name in enumerate(versions):
        r = results[name]
        n_total = r["t_stats"].shape[0]
        leaked = r["univariate"]
        pct = leaked / n_total * 100 if n_total > 0 else 0
        yval = leaked
        ax.bar(name, yval, color=colors[i], width=0.5, edgecolor="white", linewidth=0.5)
        ax.text(i, yval / 2,
                "{:,} of {:,}\n({:.1f}%)\n-- {} --".format(
                    int(leaked), n_total, pct, "LEAKING" if leaked > 0 else "SAFE"),
                ha="center", va="center", fontsize=8, color="white" if leaked > 0 else "black", fontweight="bold")
    ax.set_ylim(0, None)
    if args.order == 2:
        ax.set_ylabel("Leaked trace features")
        ax.set_title("How many trace features leaked? (simple univariate test)")
    else:
        ax.set_ylabel("Leaked trace points")
        ax.set_title("How many individual trace points leaked? (simple test)")

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
    if args.order == 2:
        ax.set_title("How many trace windows leaked? (advanced multivariate test)")
    else:
        ax.set_title("How many trace windows leaked? (advanced test)")

    if args.order == 2:
        fig.suptitle(
            "Kyber -- Summary for the Common Layman\n"
            "Can the masked version protect your secrets from second-order side-channel attacks?",
            fontsize=11, fontweight="bold", y=0.96)
    else:
        fig.suptitle(
            "Kyber -- Summary for the Common Layman\n"
            "Can the masked version protect your secrets from side-channel attacks?",
            fontsize=11, fontweight="bold", y=0.96)

    ax.set_xticks(np.arange(len(versions)))

    fig.tight_layout(rect=[0, 0, 0.97, 0.93])
    return fig


def _make_table_figure(results, args):
    """Save a compact table summarising leakage detection results.

    4 rows: unmasked/masked × univariate/multivariate.
    Columns: Implementation, Method, Statistical Test, Threshold, Detected Leakage.
    """
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties

    threshold = tvla.DEFAULT_TVLA_THRESHOLD
    alpha = tvla.DEFAULT_ALPHA
    order2 = args.order == 2

    rows_data = []
    for name, r in results.items():
        # Univariate row
        uni_label = "Welch's t-test"
        uni_thresh = "|t| > {:.2f}".format(threshold)
        uni_leak = "{:,} features".format(r["univariate"])
        method = "Univariate"
        rows_data.append((name, method, uni_label, uni_thresh, uni_leak))

        # Multivariate row
        mv_label = "Hotelling's T\u00b2"
        nu = 2 * args.nb_fixed - (args.window_size if not order2 else (args.second_order_window_size * (args.second_order_window_size + 1) // 2)) - 1
        nu = max(nu, 1)
        d = args.second_order_window_size * (args.second_order_window_size + 1) // 2 if order2 else args.window_size
        alpha_str = "{:.3e}".format(alpha)
        mv_thresh = "\u03b1 = {}".format(alpha_str)
        mv_leak = "{:,} windows (of {:,})".format(r["multivariate"], r["n_windows"])
        method = "Multivariate"
        rows_data.append((name, method, mv_label, mv_thresh, mv_leak))

    col_labels = ["Implementation", "Method", "Statistical Test", "Threshold", "Detected Leakage"]
    num_rows = len(rows_data) + 1
    row_colors = ["#f5f5f5", "#ffffff"] * 10

    grid = [col_labels] + [[cell for cell in row] for row in rows_data]

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis("off")

    col_widths = [0.14, 0.14, 0.22, 0.18, 0.32]
    table = ax.table(cellText=grid, cellLoc="center", colWidths=col_widths, loc="center")
    table.auto_set_font_size(False)

    cell_h = 0.16
    bold_font = FontProperties(weight="bold")
    normal_font = FontProperties()
    for (r, c), cell in table.get_celld().items():
        cell.set_height(cell_h)
        cell.set_text_props(fontproperties=normal_font, fontsize=9)
        if r == 0:
            cell.set_text_props(fontproperties=bold_font, fontsize=9.5, color="white")
            cell.set_facecolor("#2c3e50")
            cell.set_height(cell_h * 1.15)
        else:
            cell.set_facecolor(row_colors[r - 1])
            if c == 0:
                cell.set_text_props(fontweight="bold")
            if c == 4:
                cell.set_text_props(fontweight="bold")
        cell.set_edgecolor("#cccccc")
        cell.set_linewidth(0.5)

    fig.subplots_adjust(left=0.04, right=0.96, top=0.88, bottom=0.1)
    return fig


def make_plots(results, args):
    """Save the technical detail plot (|t| per point/feature + T² per window)."""
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
        # --- |t| per point/feature ---
        ax_t = axes[0, col]
        ax_t.plot(np.arange(r["t_stats"].shape[0]), np.abs(r["t_stats"]),
                  color="C0", lw=0.5)
        ax_t.axhline(threshold, color="r", ls="--", lw=1,
                     label="Threshold = {}".format(threshold))
        if args.order == 2:
            ax_t.set_title("{}: {} leaking features".format(name, r["univariate"]))
            ax_t.set_ylabel("Feature |t|")
        else:
            ax_t.set_title("{}: {} |t|>{} points".format(name, r["univariate"], threshold))
            ax_t.set_ylabel("|t|")
            ax_t.set_yscale("symlog")
        ax_t.legend(fontsize=7)

        # --- T² per window ---
        ax_t2 = axes[1, col]
        # colour-coded: green below threshold, red above
        t2_vals = r["t2"]

        if args.order == 2:
            # Second-order: d = features_per_window = W*(W+1)/2
            W = args.second_order_window_size
            d = W * (W + 1) // 2
        else:
            # First-order: d = window_size in points
            d = args.window_size

        n0 = args.nb_fixed
        n1 = args.nb_fixed
        nu = n0 + n1 - d - 1
        scale = (n0 + n1 - 2) * d / nu if nu > 0 else 1.0
        f_crit = tvla.f_critical(d, nu - d + 1 if nu - d + 1 > 0 else 1,
                                 tvla.DEFAULT_ALPHA)
        f_obs = t2_vals / scale

        below = f_obs < f_crit
        above = ~below
        max_t2 = max(t2_vals.max(), f_crit * scale) if t2_vals.size > 0 else f_crit * scale

        ax_t2.scatter(np.where(below)[0], t2_vals[below],
                      color="#2ecc71", s=1.5, edgecolors="none", label="Below threshold")
        ax_t2.scatter(np.where(above)[0], t2_vals[above],
                      color="#e74c3c", s=1.5, edgecolors="none", label="Above threshold")
        ax_t2.axhline(f_crit * scale, color="gray", ls="--", lw=0.8,
                      label="T² threshold = {:.1f}".format(f_crit * scale))
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

    # Save the comparison table figure
    if not args.no_table:
        table_fig = _make_table_figure(results, args)
        table_out = os.path.join(REPO_ROOT, "tvla_table.png")
        table_fig.savefig(table_out, dpi=150)
        plt.close(table_fig)
        print("Saved table to {}".format(table_out))


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Fixed-vs-random TVLA leakage study of masked vs unmasked Kyber.")
    parser.add_argument("--nb-fixed", type=int, default=256,
                        help="number of traces with the fixed (constant) secret (default 256)")
    parser.add_argument("--nb-random", type=int, default=256,
                        help="number of traces with random secrets (default 256)")
    parser.add_argument("--window-size", type=int, default=10,
                        help="multivariate window size in points (default 10)")
    parser.add_argument("--fixed-value", type=int, default=0,
                        help="constant secret value for the fixed traces (default 0)")
    parser.add_argument("--verbose", action="store_true",
                        help="print ELMO progress and -fvr output")
    parser.add_argument("--save-stats", action="store_true",
                        help="save per-version t-stats / T^2 to tvla_stats.npz")
    parser.add_argument("--no-plots", action="store_true",
                        help="do not save tvla_results.png")
    parser.add_argument("--no-table", action="store_true",
                        help="do not save tvla_table.png")
    parser.add_argument("--order", type=int, choices=[1, 2], default=1,
                        help="TVLA order: 1=first-order (default), 2=second-order centered-product")
    parser.add_argument("--second-order-window-size", type=int, default=10,
                        help="second-order centered-product window size in trace points (default 10)")
    parser.add_argument("--inject-masked-joint-leak", type=float, default=0.0,
                        metavar="SIGMA",
                        help="inject joint-only second-order leak into masked group "
                             "with per_feature_sigma=SIGMA; 0.0 = disabled (for debugging/verifying multivariate)")
    parser.add_argument("--joint-leak-fraction", type=float, default=0.1,
                        help="fraction of windows to inject when --inject-masked-joint-leak is set (default 0.1)")
    parser.add_argument("--plot-trace", action="store_true",
                        help="plot Welch univariate TVLA traces: unmasked_trace.png (unmasked), "
                             "masked_trace.png (first-order masked), masked_second_order_trace.png (second-order masked)")
    parser.add_argument("--plot-second-order-t2", action="store_true",
                        help="plot Hotelling T² multivariate TVLA traces for second-order: unmasked_second_order_t2_trace.png, masked_second_order_t2_trace.png")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.nb_fixed != args.nb_random:
        raise SystemExit(
            "--nb-fixed must equal --nb-random (balanced fixed-vs-random TVLA).")

    challenges = make_challenges(args.nb_fixed, args.nb_random, args.fixed_value)

    results = {}
    second_order_inject_meta = {}

    for label, classname, bin_path in VERSIONS:
        print("=== {} ({}) [order={}] ===".format(label, classname, args.order))
        sim = run_simulation(classname, challenges, args.verbose)
        filenames = sim.get_results_filenames()
        print("  loading {} traces...".format(len(filenames)))
        traces = load_traces(filenames)

        group0 = traces[:args.nb_fixed]
        group1 = traces[args.nb_fixed:args.nb_fixed + args.nb_fixed]
        trace_length = traces.shape[1]
        del traces

        if args.order == 1:
            # -------- FIRST-ORDER PIPELINE --------
            print("  running ELMO -fvr (ground truth)...")
            ground_truth = run_elmo_fvr(bin_path, args.nb_fixed, args.verbose)

            univariate, t_stats = tvla.count_univariate_leakages(group0, group1)
            multivariate, t2, n_singular, f_crit = tvla.count_multivariate_leakages(
                group0, group1, args.window_size)

            results[label] = {
                "ground_truth": ground_truth,
                "univariate": univariate,
                "multivariate": multivariate,
                "t_stats": t_stats,
                "t2": t2,
                "n_singular": n_singular,
                "n_windows": t2.shape[0],
                "order": 1,
            }
            print("  ground truth: {} | univariate: {} | multivariate: {}"
                  .format(ground_truth, univariate, multivariate))

            # Plot Welch unmasked univariate TVLA trace
            if args.plot_trace and label == "unmasked":
                trace_plot.plot_welch_trace(
                    t_stats,
                    group0[0],
                    os.path.join(REPO_ROOT, "unmasked_trace.png"),
                    title_prefix="Welch Unmasked Univariate TVLA",
                    y_label="Welch t-statistic",
                    show_raw_trace=True,
                    roi_length=3000)

            # Plot Welch masked univariate TVLA trace (first-order)
            if args.plot_trace and label == "masked":
                trace_plot.plot_welch_trace(
                    t_stats,
                    group0[0],
                    os.path.join(REPO_ROOT, "masked_trace.png"),
                    title_prefix="Welch Masked Univariate TVLA (First-Order)",
                    y_label="Welch t-statistic",
                    show_raw_trace=True,
                    roi_length=3000)

        else:
            # -------- SECOND-ORDER PIPELINE --------
            print("  pre-processing second-order features...")
            W = args.second_order_window_size
            features_per_window = W * (W + 1) // 2
            n_windows = trace_length // W
            total_features = n_windows * features_per_window
            print("    {} features (window_size={}, traces per group=256)".format(
                total_features, W))

            f0, f1 = tvla.preprocess_second_order_traces(
                group0, group1, args.second_order_window_size)

            # OPTIONAL: inject joint-only leak into masked group
            joint_leak_meta = {}
            if args.inject_masked_joint_leak > 0.0 and label == "masked":
                n_total_windows = f1.shape[1] // features_per_window
                f1, inject_indices = tvla.inject_second_order_joint_leak_with_indices(
                    group0, f1, args.second_order_window_size,
                    per_feature_sigma=args.inject_masked_joint_leak,
                    fraction=args.joint_leak_fraction)
                joint_leak_meta = {
                    "sigma": args.inject_masked_joint_leak,
                    "fraction": args.joint_leak_fraction,
                    "inject_indices": inject_indices,
                    "n_injected": len(inject_indices),
                    "features_per_window": features_per_window,
                    "n_total_windows": n_total_windows,
                }
                print("  INJECTED JOINT-ONLY second-order leak into MASKED: "
                      "{}/{} windows (sigma={:.2f}, fraction={:.2f})".format(
                          joint_leak_meta["n_injected"], n_total_windows,
                          args.inject_masked_joint_leak, args.joint_leak_fraction))
                second_order_inject_meta[label] = joint_leak_meta

            print("  computing second-order univariate leakages...")
            univariate, t_stats = tvla.count_second_order_univariate_leakages(
                f0, f1, args.second_order_window_size)
            print("    {} second-order features flagged".format(univariate))

            print("  computing second-order multivariate leakages...")
            multivariate, t2, n_singular, f_crit_val = tvla.count_second_order_multivariate_leakages(
                f0, f1, args.second_order_window_size)

            # Store f_crit for audit use
            if label == "masked" and joint_leak_meta:
                joint_leak_meta["f_crit"] = f_crit_val

            results[label] = {
                "univariate": univariate,
                "multivariate": multivariate,
                "t_stats": t_stats,        # per-feature Welch |t|
                "t2": t2,                  # per-window T²
                "n_singular": n_singular,
                "n_windows": t2.shape[0],
                "n_features": t_stats.shape[0],
                "f_crit": f_crit_val,
                "order": 2,
            }

            # Plot Welch unmasked second-order univariate TVLA trace
            if args.plot_trace and label == "unmasked":
                trace_plot.plot_welch_trace(
                    t_stats,
                    raw_trace=None,
                    output_path=os.path.join(REPO_ROOT, "unmasked_second_order_trace.png"),
                    title_prefix="Welch Unmasked Second-Order Univariate TVLA",
                    y_label="Feature |t|",
                    show_raw_trace=False,
                    roi_length=15000)

            # Plot Welch masked second-order univariate TVLA trace
            if args.plot_trace and label == "masked":
                trace_plot.plot_welch_trace(
                    t_stats,
                    raw_trace=None,
                    output_path=os.path.join(REPO_ROOT, "masked_second_order_trace.png"),
                    title_prefix="Welch Masked Second-Order Univariate TVLA",
                    y_label="Feature |t|",
                    show_raw_trace=False,
                    roi_length=15000)

            # Plot Hotelling T² second-order multivariate TVLA trace
            if args.plot_second_order_t2:
                W = args.second_order_window_size
                d = W * (W + 1) // 2
                n0 = args.nb_fixed
                n1 = args.nb_fixed
                nu = n0 + n1 - d - 1
                if nu > 0:
                    scale = (n0 + n1 - 2) * d / nu
                else:
                    scale = 1.0
                f_crit = tvla.f_critical(d, nu - d + 1 if nu - d + 1 > 0 else 1,
                                         tvla.DEFAULT_ALPHA)
                t2_threshold = f_crit * scale

                t_stats_label = "Welch |t|"

                if label == "unmasked":
                    output_path = os.path.join(REPO_ROOT, "unmasked_second_order_t2_trace.png")
                    t_label = "Hotelling T² Unmasked Second-Order TVLA"
                else:  # masked
                    output_path = os.path.join(REPO_ROOT, "masked_second_order_t2_trace.png")
                    t_label = "Hotelling T² Masked Second-Order TVLA"

                trace_plot.plot_hotelling_t2_trace(
                    t2=t2,
                    t2_threshold=t2_threshold,
                    t2_label=t_label,
                    t_stats=t_stats,
                    t_stats_threshold=tvla.DEFAULT_TVLA_THRESHOLD,
                    t_stats_label=t_stats_label,
                    output_path=output_path,
                    roi_length=5000)

    print_table(results, args)

    # Print joint-leak audit if injection was requested
    if second_order_inject_meta.get("masked"):
        print_joint_leak_audit(results["masked"], second_order_inject_meta["masked"], args)

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
