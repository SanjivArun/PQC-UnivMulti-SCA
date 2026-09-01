#!/usr/bin/env python3
"""Regenerate tvla_results.png from saved tvla_stats.npz (no ELMO run needed)."""

import os
import numpy as np
import tvla

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib not available")
    exit(1)

# Prevent Agg overflow ("Exceeded cell block limit") when drawing very long
# line paths under symlog scales with a clipped y-range.
matplotlib.rcParams["agg.path.chunksize"] = 20000

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
NPZ_PATH = os.path.join(REPO_ROOT, "tvla_stats.npz")

# Load saved stats
data = np.load(NPZ_PATH)
unmasked_t = data["unmasked_t"]
masked_t = data["masked_t"]
unmasked_t2 = data["unmasked_t2"]
masked_t2 = data["masked_t2"]

threshold = tvla.DEFAULT_TVLA_THRESHOLD
alpha = tvla.DEFAULT_ALPHA

# Reconstruct results dict like make_plots expects
results = {
    "unmasked": {
        "t_stats": unmasked_t,
        "t2": unmasked_t2,
        "univariate": int(np.count_nonzero(np.abs(unmasked_t) > threshold)),
        "multivariate": None,  # will compute below
        "n_windows": unmasked_t2.shape[0],
    },
    "masked": {
        "t_stats": masked_t,
        "t2": masked_t2,
        "univariate": int(np.count_nonzero(np.abs(masked_t) > threshold)),
        "multivariate": None,
        "n_windows": masked_t2.shape[0],
    }
}

# Compute multivariate leak counts (need n0, n1, window_size/d)
# This .npz comes from a SECOND-ORDER run (--order 2, --second-order-window-size 10):
#   55 centered-product features per window (W*(W+1)/2).
# Each t2 entry is one per-trace-window Hotelling T^2 over d = 55 features.
n0 = n1 = 256            # --nb-fixed default
W = 10                   # --second-order-window-size default
d = W * (W + 1) // 2     # 55 features per window
nu = n0 + n1 - d - 1     # 456 denominator degrees of freedom
scale = (n0 + n1 - 2) * d / nu
f_crit = tvla.f_critical(d, nu, alpha)

for name in ["unmasked", "masked"]:
    t2 = results[name]["t2"]
    f_obs = t2 / scale
    results[name]["multivariate"] = int(np.count_nonzero(f_obs > f_crit))

# ---- Plotting: 4 separate images ----
for name, r in results.items():
    # --- 1. Univariate: |t| per point ---
    fig_t, ax_t = plt.subplots(figsize=(10, 4))
    ax_t.plot(np.arange(r["t_stats"].shape[0]), np.abs(r["t_stats"]),
              color="C0", lw=0.5)
    ax_t.axhline(threshold, color="r", ls="--", lw=1,
                 label="Threshold = {}".format(threshold))
    ax_t.set_title("{}: Welch's t-Test Second-Order Testing ({} leaking)".format(
        name.capitalize(), r["univariate"]))
    ax_t.set_xlabel("Feature index")
    ax_t.set_ylabel("|t|")
    ax_t.set_yscale("symlog")
    ax_t.set_ylim(bottom=0, top=100)  # zoom to relevant range so data isn't compressed
    ax_t.set_xlim(left=0, right=r["t_stats"].shape[0] - 1)            # start x-axis at point 0
    ax_t.legend(fontsize=9)
    fig_t.tight_layout()
    out_t = os.path.join(REPO_ROOT, "tvla_{}_univariate.png".format(name))
    fig_t.savefig(out_t, dpi=150)
    plt.close(fig_t)
    print("Saved {}".format(out_t))

    # --- 2. Multivariate: T^2 per window ---
    fig_t2, ax_t2 = plt.subplots(figsize=(10, 4))
    t2_vals = r["t2"]
    f_obs = t2_vals / scale
    below = f_obs < f_crit
    above = ~below
    max_t2 = max(t2_vals.max(), f_crit * scale) if t2_vals.size > 0 else f_crit * scale

    ax_t2.scatter(np.where(below)[0], t2_vals[below],
                  color="#2ecc71", s=1.5, edgecolors="none", label="Below threshold")
    ax_t2.scatter(np.where(above)[0], t2_vals[above],
                  color="#e74c3c", s=1.5, edgecolors="none", label="Above threshold")
    ax_t2.axhline(f_crit * scale, color="gray", ls="--", lw=0.8,
                  label="T^2 threshold = {:.1e}".format(f_crit * scale))
    ax_t2.set_title("{}: Second-order Hotelling T\u00b2 per window ({} leaking)".format(
        name.capitalize(), r["multivariate"]))
    ax_t2.set_xlabel("Window index")
    ax_t2.set_ylabel("T\u00b2")
    ax_t2.set_yscale("symlog", linthresh=max(f_crit * scale * 1e-9, 1.0))
    ax_t2.set_ylim(bottom=0)
    ax_t2.legend(fontsize=9)
    fig_t2.tight_layout()
    out_t2 = os.path.join(REPO_ROOT, "tvla_{}_multivariate.png".format(name))
    fig_t2.savefig(out_t2, dpi=150)
    plt.close(fig_t2)
    print("Saved {}".format(out_t2))

print("\nDone: 4 images saved (second-order TVLA, W=10, 55 features/window):")
for name in ["unmasked", "masked"]:
    print("  tvla_{}_univariate.png".format(name))
    print("  tvla_{}_multivariate.png".format(name))