#!/usr/bin/env python3
"""Regenerate the TVLA table PNGs (tvla_table.png and tvla_table_simple.png)
from a saved tvla_stats.npz -- no ELMO / simulation run needed.

The .npz must come from a *second-order* run (--order 2,
--second-order-window-size 10, --nb-fixed 256 --nb-random 256) produced by
tvla_analysis.py.  That is the layout the saved stats are written in
(see `main()` in tvla_analysis.py and plot_from_npz.py):

    unmasked_t / masked_t     per-feature Welch |t|  (n_windows * 55)
    unmasked_t2 / masked_t2   per-window Hotelling T^2

Usage:

    .venv/bin/python3 table_from_npz.py
"""

import os

import numpy as np

import tvla
import tvla_analysis as ta

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
NPZ_PATH = os.path.join(REPO_ROOT, "tvla_stats.npz")

# --- Parameters matching the npz provenance (tvla_analysis defaults) ---------
N = 256            # --nb-fixed / --nb-random
W = 10             # --second-order-window-size
ORDER = 2          # second-order run

threshold = tvla.DEFAULT_TVLA_THRESHOLD
alpha = tvla.DEFAULT_ALPHA
d = W * (W + 1) // 2              # 55 centered-product features per window
nu = N + N - d - 1                # denominator degrees of freedom = 456
scale = (N + N - 2) * d / nu      # T^2 -> F scaling factor
f_crit = tvla.f_critical(d, nu, alpha)


class Args:
    nb_fixed = N
    nb_random = N
    window_size = W
    second_order_window_size = W
    order = ORDER


def build_results(data):
    """Reconstruct the per-version results dict the table functions expect."""
    results = {}
    for name in ("unmasked", "masked"):
        t_stats = data["{}_t".format(name)]
        t2 = data["{}_t2".format(name)]
        n_windows = t2.shape[0]
        f_obs = t2 / scale
        results[name] = {
            "univariate": int(np.count_nonzero(np.abs(t_stats) > threshold)),
            "multivariate": int(np.count_nonzero(f_obs > f_crit)),
            "t_stats": t_stats,
            "t2": t2,
            "n_singular": 0,
            "n_windows": n_windows,
            "trace_length": n_windows * W,
            "order": ORDER,
        }
    return results


def main():
    data = np.load(NPZ_PATH)
    results = build_results(data)
    data.close()

    for name, r in results.items():
        print("  {}: univariate={:,} multivariate={:,} windows={:,} "
              "trace_length={:,}".format(
                  name, r["univariate"], r["multivariate"], r["n_windows"],
                  r["trace_length"]))

    args = Args()

    table_fig = ta._make_table_figure(results, args)
    table_out = os.path.join(REPO_ROOT, "tvla_table.png")
    table_fig.savefig(table_out, dpi=150)
    print("Saved {}".format(table_out))

    simple_fig = ta._make_simple_table_figure(results, args)
    simple_out = os.path.join(REPO_ROOT, "tvla_table_simple.png")
    simple_fig.savefig(simple_out, dpi=150)
    print("Saved {}".format(simple_out))

    print("\nDone: tables regenerated from {}".format(NPZ_PATH))


if __name__ == "__main__":
    main()
