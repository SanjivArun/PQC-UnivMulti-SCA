#!/usr/bin/env python3
"""trace_plot.py -- Plot Welch Unmasked Univariate TVLA trace data from ELMO.

Called by tvla_analysis.py via:   trace_plot.plot_unmasked_welch_trace(...)

Produces:  unmasked_trace.png (Welch t-stat line graph + sample power trace overlay)
"""


def plot_unmasked_welch_trace(t_stats, raw_trace, output_path, threshold=4.5):
    """Save the Welch t-statistic as a line-graph PNG with a sample power trace overlay.

    Parameters
    ----------
    t_stats : np.ndarray          - per-point Welch t-statistic, shape (n_points,)
    raw_trace : np.ndarray        - one raw power trace from the ELMO output, shape (n_points,)
    output_path : str             - PNG output file path
    threshold : float             - |t| threshold (default 4.5)
    """
    import numpy as np

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping trace plot.")
        return

    n_points = len(t_stats)
    positions = np.arange(n_points)

    fig, (ax_t, ax_raw) = plt.subplots(
        2, 1, figsize=(14, 5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )
    fig.suptitle("Welch Unmasked Univariate TVLA", fontsize=13, fontweight="bold")

    # --- Panel A: Welch t-statistic line graph ---
    ax_t.plot(positions, t_stats, color="#1f77b4", lw=0.8, zorder=5, label="Welch t")

    ax_t.axhline(0, color="silver", ls="-", lw=1, zorder=1)
    ax_t.axhline(
        threshold,
        color="red", ls="--", lw=1.2,
        label="|t| > {:.1f}".format(threshold),
        zorder=2,
    )
    ax_t.axhline(-threshold, color="red", ls="--", lw=1.2, zorder=2)
    ax_t.axhspan(-threshold, threshold, color="green", alpha=0.06, zorder=0)
    ax_t.set_ylabel("Welch t-statistic")
    ax_t.legend(loc="upper right", fontsize=9)
    ax_t.set_ylim(-7, 7)

    # --- Panel B: Sample raw power trace ---
    ax_raw.fill_between(positions, raw_trace, alpha=0.3, color="gray")
    ax_raw.set_xlabel("Trace position")
    ax_raw.set_ylabel("Power")

    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved trace plot to {output_path}")
