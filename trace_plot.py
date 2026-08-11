#!/usr/bin/env python3
"""trace_plot.py -- Plot Welch Univariate TVLA trace data from ELMO.

Called by tvla_analysis.py via:   trace_plot.plot_unmasked_welch_trace(...)
                                  trace_plot.plot_welch_trace(...)

Produces:  unmasked_trace.png, masked_trace.png (Welch t-stat line graph +
           sample power trace overlay).
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


def plot_welch_trace(t_stats, raw_trace, output_path, threshold=4.5,
                     title_prefix="Welch Masked Univariate TVLA",
                     y_label="Welch t-statistic",
                     show_raw_trace=True,
                     roi_length=None,
                     normalize_trace=False):
    """Save the Welch t-statistic as a line-graph PNG with a sample power trace overlay.

    Parameters
    ----------
    t_stats : np.ndarray          - per-point or per-feature Welch t-statistic, shape (n,)
    raw_trace : np.ndarray        - one raw power trace from the ELMO output, shape (n,)
    output_path : str             - PNG output file path
    threshold : float             - |t| threshold (default 4.5)
    title_prefix : str            - prefix for the plot title
    y_label : str                 - label for the y-axis of the t-statistic panel
    show_raw_trace : bool         - whether to show Panel B (raw power trace).
                                    Set False for second-order plots (where raw
                                    trace overlay doesn't apply).
    roi_length : int or None      - number of leading points to display.  If
                                    *None* (default) the full trace is shown.
                                    Setting to ~5000 zooms into the region of
                                    interest for readability.
    normalize_trace : bool        - mean-centre and divide by std-dev so the
                                    power y-axis stays in a readable range
                                    (e.g. ±5 for the masked reference).
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

    # ROI slice: limit x-axis to the first *roi_length* points so that the
    # readable region is comparable to the masked reference plot.
    end_idx = min(n_points, roi_length) if roi_length is not None else n_points

    pos_display = np.arange(end_idx, dtype=np.int64)
    t_display = t_stats[:end_idx]
    if raw_trace is not None:
        raw_display = raw_trace[:end_idx]
        if normalize_trace:
            raw_display = (raw_display - raw_display.mean()) / raw_display.std()
    else:
        raw_display = None

    if show_raw_trace:
        fig, (ax_t, ax_raw) = plt.subplots(
            2, 1, figsize=(14, 5), sharex=True,
            gridspec_kw={"height_ratios": [3, 1]}
        )
    else:
        fig, ax_t = plt.subplots(1, 1, figsize=(14, 4.5))
        ax_raw = None

    fig.suptitle(title_prefix, fontsize=13, fontweight="bold")

    # --- Panel A: Welch t-statistic line graph ---
    # Use alpha transparency so overlapping points show density
    ax_t.plot(pos_display, t_display, color="#1f77b4", lw=0.8, zorder=5, label="Welch t", alpha=0.9)

    ax_t.axhline(0, color="silver", ls="-", lw=1, zorder=1)
    ax_t.axhline(
        threshold,
        color="red", ls="--", lw=1.2,
        label="|t| > {:.1f}".format(threshold),
        zorder=2,
    )
    ax_t.axhline(-threshold, color="red", ls="--", lw=1.2, zorder=2)
    ax_t.axhspan(-threshold, threshold, color="green", alpha=0.06, zorder=0)
    ax_t.set_ylabel(y_label)
    max_t = float(np.nanmax(np.abs(t_display)))
    if np.isinf(max_t) or np.isnan(max_t):
        max_t = threshold
    y_max_display = max(threshold, max_t)
    ax_t.set_ylim(-y_max_display * 1.12, y_max_display * 1.12)

    # --- Panel B: Sample raw power trace ---
    if show_raw_trace and raw_display is not None:
        ax_raw.plot(pos_display, raw_display, color="gray", lw=0.3, zorder=1, alpha=0.4)
        ax_raw.set_xlabel("Trace position")
        ax_raw.set_ylabel("Power")

    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved trace plot to {output_path}")
