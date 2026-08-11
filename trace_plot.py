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


def plot_hotelling_t2_trace(t2, t2_threshold, t2_label=None,
                           t_stats=None, t_stats_threshold=4.5,
                           t_stats_label=None,
                           output_path="/tmp/hotelling_t2.png",
                           roi_length=None):
    """Save Hotelling T² multivariate TVLA trace data as a line-graph PNG.

    Two-panel figure:
      Panel A: T² per window with threshold line and green "safe" zone
      Panel B: Welch |t| per-feature (second-order), if provided

    Parameters
    ----------
    t2 : np.ndarray
        Per-window T² values, shape (n_windows,).
    t2_threshold : float
        F critical value scaled to T² units (t2_threshold = f_crit * scale).
        Windows above this line are above the detection threshold.
    t2_label : str or None
        Label for the T² panel y-axis (defaults to "T²").
    t_stats : np.ndarray or None
        Per-feature Welch absolute t-statistics, shape (n_features,).
        If provided, Panel B shows |t| for comparison.
    t_stats_threshold : float
        |t| threshold, default 4.5.
    t_stats_label : str or None
        Label for the Welch |t| panel y-axis (defaults to "Welch |t|").
    output_path : str
        PNG output file path.
    roi_length : int or None
        Number of leading windows to display in Panel A. If *None*,
        the full result is shown.
    """
    import numpy as np

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping trace plot.")
        return

    n_windows = len(t2)
    pos = np.arange(n_windows)

    # ROI slicing
    end = min(n_windows, roi_length) if roi_length is not None else n_windows
    pos_display = pos[:end]
    t2_display = t2[:end]

    # --- Panel A: T² line graph ---
    fig, ax_t2 = plt.subplots(1, 1, figsize=(14, 3.5))
    fig.suptitle("Hotelling T² Multivariate TVLA", fontsize=13, fontweight="bold")

    ax_t2.plot(pos_display, t2_display, color="#1f77b4", lw=0.8, zorder=5, label="T²")

    # Green zone below threshold
    ax_t2.axhspan(0, t2_threshold, color="green", alpha=0.08, zorder=0)
    # Red threshold line
    ax_t2.axhline(t2_threshold, color="red", ls="--", lw=1.2,
                  label="T² threshold = {:.1f}".format(t2_threshold), zorder=2)

    # Mark leaking windows (above threshold)
    leaking_mask = t2_display > t2_threshold
    if leaking_mask.any():
        ax_t2.scatter(pos_display[leaking_mask], t2_display[leaking_mask],
                      color="red", s=8, zorder=10, label="Leaking windows",
                      edgecolors="none")

    ax_t2.set_ylabel(t2_label if t2_label else "T²")
    ax_t2.set_xlabel("Window index")
    ax_t2.legend(loc="upper right", fontsize=9)

    max_t2 = max(float(t2_display.max()), t2_threshold) if t2_display.size > 0 else t2_threshold
    ax_t2.set_ylim(0, max_t2 * 1.12)

    # --- Panel B: Welch |t| if provided ---
    if t_stats is not None:
        ax_t2 = None

    if t_stats is not None:
        n_features = len(t_stats)
        end_b = min(n_features, roi_length) if roi_length is not None else n_features
        pos_display_b = np.arange(end_b)
        t_stats_display = t_stats[:end_b]

        fig, (ax_t2, ax_welch) = plt.subplots(
            2, 1, figsize=(14, 6), sharex=True,
            gridspec_kw={"height_ratios": [3, 1]}
        )
        fig.suptitle("Hotelling T² + Welch |t| Multivariate TVLA", fontsize=13, fontweight="bold")

        # Panel A: T² line graph (same as above)
        ax_t2.plot(pos_display, t2_display, color="#1f77b4", lw=0.8, zorder=5, label="T²")
        ax_t2.axhspan(0, t2_threshold, color="green", alpha=0.08, zorder=0)
        ax_t2.axhline(t2_threshold, color="red", ls="--", lw=1.2,
                      label="T² threshold = {:.1f}".format(t2_threshold), zorder=2)
        if leaking_mask.any():
            ax_t2.scatter(pos_display[leaking_mask], t2_display[leaking_mask],
                          color="red", s=8, zorder=10, label="Leaking windows",
                          edgecolors="none")
        ax_t2.set_ylabel(t2_label if t2_label else "T²")
        ax_t2.set_xlabel("Window index")
        ax_t2.legend(loc="upper right", fontsize=9)
        ax_t2.set_ylim(0, max_t2 * 1.12)

        # Panel B: Welch |t| per feature
        ax_welch.plot(pos_display_b, t_stats_display, color="#1f77b4", lw=0.8,
                      zorder=5, label="Welch |t|", alpha=0.9)
        ax_welch.axhline(t_stats_threshold, color="red", ls="--", lw=1.2,
                        label="|t| > {:.1f}".format(t_stats_threshold), zorder=2)
        ax_welch.axhline(-t_stats_threshold, color="red", ls="--", lw=1.2, zorder=2)
        ax_welch.axhspan(-t_stats_threshold, t_stats_threshold,
                        color="green", alpha=0.06, zorder=0)
        max_welch = float(np.nanmax(np.abs(t_stats_display)))
        if np.isinf(max_welch) or np.isnan(max_welch):
            max_welch = t_stats_threshold
        y_max_welch = max(t_stats_threshold, max_welch)
        ax_welch.set_ylim(-y_max_welch * 1.12, y_max_welch * 1.12)
        ax_welch.set_ylabel(t_stats_label if t_stats_label else "Welch |t|")
        ax_welch.set_xlabel("Feature index")
        ax_welch.legend(loc="upper right", fontsize=9)

    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved T² trace plot to {output_path}")
