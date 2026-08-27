"""tvla.py -- Leakage-detection statistics for ELMO power traces (pure NumPy).

This module implements the two leakage tests used in the side-channel literature
to detect *where* (and *how many*) instruction points of a power trace depend on
secret data, in the "fixed vs random" TVLA (Test Vector Leakage Assessment)
experimental design:

* **Univariate Welch's t-test** (`welch_t_statistic`, `count_univariate_leakages`):
  compares the two trace groups point by point. A point "leaks" when the absolute
  value of its t-statistic exceeds the classic TVLA threshold ``|t| > 4.5``.

* **Multivariate Hotelling's T^2 test** (`hotelling_t2_per_window`,
  `count_multivariate_leakages`): groups the trace into consecutive non-overlapping
  windows of ``window_size`` points and tests each window *jointly*. A window
  "leaks" when its T^2 exceeds the critical value of the F distribution (the exact
  finite-sample distribution of T^2 under the null hypothesis).

The two groups are matrices of shape ``(n_traces, n_points)``; both groups must
have the **same** number of points.

The module has **no dependency on ELMO** and no ``scipy`` dependency: the
F-distribution quantile is computed exactly from the regularized incomplete beta
function (Numerical-Recipes-style continued fraction), so everything runs on
plain NumPy.

Why fixed-vs-random:
    group0 = the first  N traces, produced with the *same* (fixed) input each time;
    group1 = the next  N traces, produced with *random* inputs.
    Any point whose power differs between the two groups "leaks" (depends on data).
"""

import math

import numpy as np

# --------------------------------------------------------------------------
# Thresholds / significance levels
# --------------------------------------------------------------------------

#: Standard TVLA threshold: a point leaks when |t| > 4.5 (see e.g. the
#: original TVLA methodology and ELMO's own FIXEDVSRANDOMFAIL macro).
DEFAULT_TVLA_THRESHOLD = 4.5

#: Two-sided normal significance level implied by ``|t| > DEFAULT_TVLA_THRESHOLD``.
#: Used as the default per-window alpha for the multivariate test so that the
#: two tests have matched per-test false-positive rates.
DEFAULT_ALPHA = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(DEFAULT_TVLA_THRESHOLD / math.sqrt(2.0))))


def welch_t_statistic(group0, group1):
    """Per-point Welch's t-statistic between two trace groups.

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n0, p) / (n1, p)
        Power traces of the two groups (traces on axis 0, points on axis 1).

    Returns
    -------
    np.ndarray, shape (p,)
        t = (mean0 - mean1) / sqrt(var0 / n0 + var1 / n1)  with unbiased sample
        variances (``ddof=1``). Points where the denominator is zero get:
        ``0`` if the two means coincide, otherwise ``+/- inf`` (the sign of the
        mean difference) -- mirroring ELMO's own ``getttest`` behavior.

    Notes
    -----
    In the unmasked fixed-vs-random experiment the "fixed" group is a repeated
    identical trace, so its variance is exactly 0. This is handled naturally: the
    t-statistic reduces to ``(mean0 - mean1) / sqrt(var1 / n1)``.
    """
    a = np.asarray(group0, dtype=np.float64)
    b = np.asarray(group1, dtype=np.float64)
    n0, n1 = a.shape[0], b.shape[0]
    mean0, mean1 = a.mean(axis=0), b.mean(axis=0)
    var0, var1 = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (mean0 - mean1) / np.sqrt(var0 / n0 + var1 / n1)
    # 0/0 -> nan  -> 0 (no detectable difference);  x/0 -> +/-inf (kept as inf).
    return np.nan_to_num(t, nan=0.0, posinf=np.inf, neginf=-np.inf)


def count_univariate_leakages(group0, group1, threshold=DEFAULT_TVLA_THRESHOLD):
    """Count leaking points using the univariate Welch's t-test.

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n, p)
        The two trace groups.
    threshold : float
        A point leaks when ``abs(t) > threshold`` (default 4.5).

    Returns
    -------
    (n_leakages, t_stats)
        ``n_leakages`` : int -- number of leaking points.
        ``t_stats`` : np.ndarray, shape (p,) -- the per-point t-statistics
        (kept so callers can plot / save them for debugging).
    """
    t = welch_t_statistic(group0, group1)
    return int(np.count_nonzero(np.abs(t) > threshold)), t


def hotelling_t2_per_window(group0, group1, window_size):
    """Two-sample Hotelling's T^2 for consecutive non-overlapping windows of points.

    For each window of ``window_size`` consecutive points:

        T^2 = (n0 * n1 / (n0 + n1)) * (x0bar - x1bar)^T  S_pooled^{-1} (x0bar - x1bar)

    with pooled covariance
        S_pooled = ((n0 - 1) * S0 + (n1 - 1) * S1) / (n0 + n1 - 2).

    A (near-)singular pooled covariance (e.g. constant power columns) is inverted
    with the Moore-Penrose pseudo-inverse; the number of windows where the
    covariance was rank-deficient is reported so the caller can judge how much of
    the result relies on the pseudo-inverse.

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n0, p) / (n1, p)
        The two trace groups (same number of points ``p``).
    window_size : int
        Number of consecutive points per window. Must satisfy
        ``n0 + n1 > window_size + 1`` so the pooled covariance is estimable.

    Returns
    -------
    (t2, n_singular)
        ``t2`` : np.ndarray, shape (nb_windows,) -- T^2 per window, where
        ``nb_windows = max(1, (p0 - window_size) // window_size + 1)``
        (trailing points that fall short of a full window are silently dropped).
        ``n_singular`` : int -- number of windows whose pooled covariance matrix
        had rank < window_size (handled with the pseudo-inverse).
    """
    a = np.asarray(group0, dtype=np.float64)
    b = np.asarray(group1, dtype=np.float64)
    n0, p0 = a.shape
    n1, p1 = b.shape
    if p0 != p1:
        raise ValueError("Both groups must have the same number of points: {} vs {}".format(p0, p1))
    if n0 + n1 <= window_size + 1:
        raise ValueError(
            "Need n0 + n1 > window_size + 1 to estimate the pooled covariance "
            "(got n0={}, n1={}, window_size={}).".format(n0, n1, window_size))
    if window_size < 1:
        raise ValueError("window_size must be >= 1, got {}.".format(window_size))

    nb_windows = max(1, (p0 - window_size) // window_size + 1)
    t2 = np.zeros(nb_windows, dtype=np.float64)
    n_singular = 0
    factor = n0 * n1 / (n0 + n1)

    for w in range(nb_windows):
        start = w * window_size
        end = start + window_size
        if end > p0:
            break
        x0 = a[:, start:end]
        x1 = b[:, start:end]
        mean0 = x0.mean(axis=0)
        mean1 = x1.mean(axis=0)
        cov0 = np.cov(x0, rowvar=False, ddof=1)
        cov1 = np.cov(x1, rowvar=False, ddof=1)
        pooled = ((n0 - 1) * cov0 + (n1 - 1) * cov1) / (n0 + n1 - 2)
        if np.linalg.matrix_rank(pooled) < window_size:
            n_singular += 1
        diff = mean0 - mean1
        t2[w] = factor * float(diff @ np.linalg.pinv(pooled) @ diff)

    return t2, n_singular


def f_critical(df1, df2, alpha):
    """Upper ``alpha`` quantile of the F distribution (no scipy).

    Uses the identity ``F_p(df1, df2) = (df1/df2) * x/(1 - x)`` where ``x`` is the
    ``p``-quantile of the Beta(df1/2, df2/2) distribution, with the Beta quantile
    found by bisection on the regularized incomplete beta function (computed with
    a Numerical-Recipes-style continued fraction).

    Parameters
    ----------
    df1, df2 : float
        Numerator / denominator degrees of freedom (both > 0).
    alpha : float
        Right-tail probability (``0 < alpha < 1``).

    Returns
    -------
    float
        The value ``F_crit`` with ``P(F > F_crit) = alpha``.
    """
    if df1 <= 0 or df2 <= 0:
        raise ValueError("Degrees of freedom must be positive: df1={}, df2={}".format(df1, df2))
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1), got {}.".format(alpha))
    # If X ~ F(df1, df2) then Y = (df1*X/df2) / (1 + df1*X/df2) ~ Beta(df1/2, df2/2).
    # Inverting: F_p = (df2/df1) * x/(1 - x) with x the Beta(df1/2, df2/2) p-quantile.
    a, b = df1 / 2.0, df2 / 2.0
    x = _beta_quantile(a, b, 1.0 - alpha)  # upper tail: I_x(a, b) = 1 - alpha
    return (df2 / df1) * x / (1.0 - x)


def count_multivariate_leakages(group0, group1, window_size=50, alpha=None):
    """Count leaking windows using Hotelling's T^2 per non-overlapping window.

    Under the null hypothesis the per-window statistic follows an F distribution:

        T^2  ~  ((n0 + n1 - 2) * d / (n0 + n1 - d - 1)) * F(d, n0 + n1 - d - 1)

    with ``d = window_size``. A window leaks when its T^2 exceeds the F critical
    value (so the test is exact for finite samples, not just asymptotically).

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n, p)
        The two trace groups.
    window_size : int
        Number of consecutive points per window (default 50).
    alpha : float, optional
        Per-window significance level. Defaults to ``DEFAULT_ALPHA`` -- the
        two-sided level matched to the univariate ``|t| > 4.5`` threshold, so both
        tests share the same per-test false-positive rate.

    Returns
    -------
    (n_leakages, t2, n_singular, f_crit)
        ``n_leakages`` : int -- number of leaking windows.
        ``t2`` : np.ndarray -- per-window T^2 values (for plotting / saving).
        ``n_singular`` : int -- windows whose covariance needed the pseudo-inverse.
        ``f_crit`` : float -- the F critical value used (for debugging).
    """
    if alpha is None:
        alpha = DEFAULT_ALPHA
    n0 = group0.shape[0]
    n1 = group1.shape[0]
    d = window_size
    nu = n0 + n1 - d - 1          # denominator degrees of freedom
    if nu < 1:
        raise ValueError(
            "n0 + n1 - window_size - 1 must be >= 1 "
            "(got n0={}, n1={}, window_size={}).".format(n0, n1, d))
    scale = (n0 + n1 - 2) * d / nu
    t2, n_singular = hotelling_t2_per_window(group0, group1, d)
    f_obs = t2 / scale
    f_crit = f_critical(d, nu, alpha)
    return int(np.count_nonzero(f_obs > f_crit)), t2, n_singular, f_crit


# --------------------------------------------------------------------------
# Second-order TVLA (centered-product transform + same statistical tests)
# --------------------------------------------------------------------------

def _extract_upper_triangular_indices(size):
    """Return (row_indices, col_indices) for the upper-triangular part of a size×size matrix.

    This includes the diagonal (i == j) so that squared terms are also captured.
    For size=3 the pairs are (0,0), (0,1), (0,2), (1,1), (1,2), (2,2) = 6 = 3*4/2.

    Parameters
    ----------
    size : int
        Dimension of the square window.

    Returns
    -------
    row_idx, col_idx : np.ndarray of int
        Index arrays of shape (k,) where k = size*(size+1)//2.
    """
    rows, cols = np.triu_indices(size)
    return rows, cols


def preprocess_second_order_traces(group0, group1, window_size=10):
    """Transform raw power traces into second-order features via centered products.

    This implements the standard centered-product second-order transform from the
    side-channel literature (e.g., Narasimhan, Balbuzhev, and Biba). The key idea is
    that a 2-share masked implementation may have no first-order leakage (each
    individual sample is independent of the secret), but pairwise products of
    centered samples from the same instruction can still leak because both samples
    depend on the same secret share.

    Preprocessing steps:

    1. Compute a **global** mean trace across ALL traces from BOTH groups combined.
       Under H0 (no leakage) both groups share the same distribution, so a single
       mean estimate is the correct model. Under H1, global centering is more
       conservative: if second-order leakage survives it, the signal is stronger.

    2. Center each trace by subtracting the global mean:
       x'(i) = x(i) - mu_global(i)

    3. Split each centered trace into non-overlapping windows of `window_size` points.

    4. For each window, compute the upper-triangular part of the outer product:
       For a 10-point window this yields 10*11/2 = 55 features.
       These features are the pairwise products (including self-products):
       x'(s)*x'(s), x'(s)*x'(s+1), ..., x'(s+W-1)*x'(s+W-1)

    Mathematical justification:
    If x(i) = f(s_share) + n(i) and x(j) = g(s_share) + n(j), then:
        (x(i) - mu(i)) * (x(j) - mu(j)) ≈ f(s_share) * g(s_share) + noise
    where the cross-term f(s_share)*g(s_share) depends on the secret share.
    Centering removes the DC offset so we don't pick up the masked values themselves,
    and the product averages out the independent noise terms.

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n0, p) / (n1, p)
        The two trace groups (traces on axis 0, points on axis 1).
    window_size : int
        Number of consecutive trace points per window for the centered-product
        transform (default 10). With W=10, each window yields W*(W+1)/2 = 55
        second-order features. The total number of features is
        (p // window_size) * window_size*(window_size+1)/2.

    Returns
    -------
    f0, f1 : np.ndarray, shape (n0, F) / (n1, F)
        Second-order feature matrices where each column is one pairwise product
        within a window. F = (p // window_size) * window_size*(window_size+1)//2.
        Columns are ordered by: (window_index, pair_index_within_window).

    Raises
    ------
    ValueError
        If groups have different numbers of points.

    Memory analysis:
        With p=93000 points, window_size=10:
        - Number of windows per trace: p // 10 = 9300
        - Features per window: 10*11//2 = 55
        - Total feature columns: 9300 * 55 = 511,500
        - For 1024 traces at np.float64: 1024 * 511500 * 8 bytes ≈ 4.0 GB
        - For 256 traces: 1.0 GB -- acceptable for research work
        - To reduce memory: decrease window_size (fewer features per window) or
          increase window_size (fewer windows) at the cost of temporal resolution.
    """
    a = np.asarray(group0, dtype=np.float32)
    b = np.asarray(group1, dtype=np.float32)
    n0, p0 = a.shape
    n1, p1 = b.shape
    if p0 != p1:
        raise ValueError("Both groups must have the same number of points: {} vs {}".format(p0, p1))
    if window_size < 2:
        raise ValueError("window_size must be >= 2 for second-order features (got {}).".format(window_size))

    # Step 1+2: Global centering over ALL traces from both groups.
    # Concatenate temporarily just to compute the mean; avoids full memory allocation.
    global_mean = np.concatenate([a, b], axis=0).mean(axis=0)  # shape: (p,)
    a_centered = a - global_mean  # shape: (n0, p)
    b_centered = b - global_mean  # shape: (n1, p)

    # Step 3+4: Window the centered traces and extract upper-triangular products.
    # We process each window index sequentially to avoid creating the full O(n*p²)
    # array in memory at once -- we build the output column-by-column.
    n_windows = p0 // window_size
    features_per_window = window_size * (window_size + 1) // 2
    f0 = np.empty((n0, n_windows * features_per_window), dtype=np.float32)
    f1 = np.empty((n1, n_windows * features_per_window), dtype=np.float32)

    row_idx, col_idx = _extract_upper_triangular_indices(window_size)

    for w in range(n_windows):
        start = w * window_size
        end = start + window_size
        # Window data: shape (n, window_size) for both groups
        w0 = a_centered[:, start:end]  # (n0, window_size)
        w1 = b_centered[:, start:end]  # (n1, window_size)
        # Upper-triangular products: for each trace, pick the (i,j) pairs
        # and compute x'(i) * x'(j) -- shape (n, features_per_window)
        f0[:, w * features_per_window:(w + 1) * features_per_window] = \
            w0[:, row_idx] * w0[:, col_idx]
        f1[:, w * features_per_window:(w + 1) * features_per_window] = \
            w1[:, row_idx] * w1[:, col_idx]

    del a_centered, b_centered
    return f0, f1


def count_second_order_univariate_leakages(group0, group1, window_size=10, threshold=DEFAULT_TVLA_THRESHOLD):
    """Count leaking second-order features using Welch's t-test on centered products.

    After the centered-product transform, each feature column corresponds to one
    pairwise product (i, j) within one trace window. This function applies the
    standard Welch's t-test pointwise to each such feature column -- the direct
    analogue of :func:`count_univariate_leakages` for second-order features.

    A feature "leaks" when the absolute value of its t-statistic exceeds the TVLA
    threshold ``|t| > threshold``.

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n0, F) / (n1, F)
        Second-order feature matrices (output of
        :func:`preprocess_second_order_traces`).
    window_size : int
        Window size used for the centered-product transform (used for metadata
        in the caller, not for computation since the features are already extracted).
    threshold : float
        A feature leaks when ``abs(t) > threshold`` (default 4.5).

    Returns
    -------
    (n_leakages, t_stats)
        ``n_leakages`` : int -- number of leaking feature columns.
        ``t_stats`` : np.ndarray, shape (F,) -- the per-feature t-statistics.
        The columns can be decoded as: window_index = i // features_per_window,
        pair_index_within_window = i % features_per_window.

    Notes
    -----
    This is the univariate second-order analogue of :func:`count_univariate_leakages`.
    Each feature is treated independently: no attempt is made to model the
    dependence between pairwise products within the same window. For that, use
    :func:`count_second_order_multivariate_leakages` instead.
    """
    t = welch_t_statistic(group0, group1)
    return int(np.count_nonzero(np.abs(t) > threshold)), t


def count_second_order_multivariate_leakages(group0, group1, window_size=10, alpha=None):
    """Count leaking second-order windows using Hotelling's T^2 on centered products.

    After the centered-product transform, each trace window yields
    ``window_size * (window_size + 1) // 2`` features (all pairwise products within
    the window). This function applies Hotelling's T^2 test **per trace-window**: each
    window of centered-product features is tested as a multivariate block.

    The rationale: if multiple pairwise products within the same trace window leak
    information about the secret share, testing them jointly (rather than one at a
    time) increases sensitivity. The pooled covariance matrix captures the
    correlation structure between pairs, and T^2 aggregates evidence across all
    features in the window.

    Under the null hypothesis the per-window statistic follows an F distribution:

        T^2  ~  ((n0 + n1 - 2) * d / (n0 + n1 - d - 1)) * F(d, n0 + n1 - d - 1)

    with ``d = window_size * (window_size + 1) // 2`` (the number of features per
    trace window). A window "leaks" when its T^2 exceeds the F critical value.

    Parameters
    ----------
    group0, group1 : np.ndarray, shape (n0, F) / (n1, F)
        Second-order feature matrices (output of
        :func:`preprocess_second_order_traces`).
    window_size : int
        The trace-point window size used for the centered-product transform. This
        determines ``d`` (the degrees of freedom): d = window_size*(window_size+1)//2.
    alpha : float, optional
        Per-window significance level. Defaults to ``DEFAULT_ALPHA`` -- the
        two-sided level matched to the univariate ``|t| > 4.5`` threshold.

    Returns
    -------
    (n_leakages, t2, n_singular, f_crit)
        ``n_leakages`` : int -- number of leaking trace windows.
        ``t2`` : np.ndarray, shape (n_trace_windows,) -- per-window T^2 values.
        ``n_singular`` : int -- windows whose pooled covariance needed the
            pseudo-inverse (rank-deficient).
        ``f_crit`` : float -- the F critical value used.

    Raises
    ------
    ValueError
        If the number of traces is too small to estimate the pooled covariance
        matrix (``n0 + n1 <= d + 1`` where d is the number of features per window).

    Notes
    -----
    This is the multivariate second-order analogue of
    :func:`count_multivariate_leakages`. Since the features are already grouped
    into trace windows, we pass them to :func:`hotelling_t2_per_window` with a
    window_size parameter equal to the number of features per second-order window
    (d = W*(W+1)/2). The function then treats each consecutive block of d features
    as one multivariate observation -- but since we already have one block per trace
    window, each "Hotelling window" contains exactly one trace window's worth of
    features, i.e. it collapses to a single per-trace-window T^2 test.

    If n0 + n1 equals the number of trace windows exactly (which is the common case:
    all features per window correspond to one trace-point window), then the per-
    window T^2 is computed over the entire trace window jointly.
    """
    if alpha is None:
        alpha = DEFAULT_ALPHA

    a = np.asarray(group0, dtype=np.float64)
    b = np.asarray(group1, dtype=np.float64)
    n0, F = a.shape
    n1, _ = b.shape
    if F != b.shape[1]:
        raise ValueError("Both groups must have the same number of features: {} vs {}".format(F, b.shape[1]))

    n_trace_windows = F // (window_size * (window_size + 1) // 2)
    d = window_size * (window_size + 1) // 2  # features per trace window

    nu = n0 + n1 - d - 1
    if nu < 1:
        raise ValueError(
            "Not enough traces for second-order multivariate test: "
            "need n0 + n1 > d + 1 (got n0={}, n1={}, features_per_window d={}*({}+1)//2={}). "
            "Reduce --second-order-window-size or increase --nb-fixed/--nb-random.".format(
                n0, n1, window_size, window_size, d))
    # trace-window. We test each trace-window jointly. Since the features are already
    # arranged as [w0_features, w1_features, ..., wN_features], we use
    # window_size=d in hotelling_t2_per_window so that ONE "Hotelling window" = the
    # entire feature vector (one feature block at a time), giving one T^2 per
    # trace-window of original trace points.
    t2, n_singular = hotelling_t2_per_window(a, b, d)
    nu = n0 + n1 - d - 1
    if nu < 1:
        raise ValueError(
            "n0 + n1 - d - 1 must be >= 1 to estimate T^2 "
            "(got n0={}, n1={}, window_size={}).".format(n0, n1, window_size))
    scale = (n0 + n1 - 2) * d / nu
    f_obs = t2 / scale
    f_crit = f_critical(d, nu, alpha)
    return int(np.count_nonzero(f_obs > f_crit)), t2, n_singular, f_crit


def inject_second_order_leak(f1, window_size, magnitude=0.5, fraction=0.1,
                             window_index=None, seed=0):
    """Inject an intentional second-order leak into the random group's features.

    For a chosen subset of trace windows, adds a constant bias to ALL
    centered-product features of that window in the random-group feature matrix
    ``f1``. The bias is expressed **relative to each feature's own standard
    deviation** -- a bias of ``magnitude`` standard deviations in every chosen
    feature -- so the same ``magnitude`` value is meaningful regardless of the
    absolute scale of the power traces (ELMO power values can be in the thousands,
    so an absolute bias of 0.5 would be invisible).

    This mimics a second-order mean bias (E[c(i)c(j)] differs between the fixed
    and the random groups) and is used to validate that the second-order
    univariate / multivariate tests actually detect leakage when it is present.

    Parameters
    ----------
    f1 : np.ndarray, shape (n1, F)
        Random-group second-order feature matrix (output of
        :func:`preprocess_second_order_traces`). Modified in place.
    window_size : int
        Trace-point window size; ``F = n_windows * window_size*(window_size+1)//2``.
    magnitude : float
        Bias in **standard deviations** added to every feature of each chosen
        window. Must be ``> 0``. ``1.0`` = a one-sigma shift, which is detectable.
    fraction : float
        Fraction of windows to leak (``0 < fraction <= 1``). Ignored when
        ``window_index`` is given.
    window_index : int, optional
        If given, leak exactly this window instead (must be in ``[0, n_windows)``).
    seed : int
        Seed for the RNG that picks which windows are leaked when using
        ``fraction``, so runs are reproducible.

    Returns
    -------
    (f1, n_injected)
        ``f1`` (the same array, mutated) and the number of windows injected with
        the leak.
    """
    features_per_window = window_size * (window_size + 1) // 2
    n_windows = f1.shape[1] // features_per_window
    if magnitude <= 0:
        raise ValueError("magnitude must be > 0, got {}.".format(magnitude))

    if window_index is None:
        if not 0 < fraction <= 1:
            raise ValueError("fraction must be in (0, 1], got {}.".format(fraction))
        n_inject = int(round(n_windows * fraction))
        rng = np.random.RandomState(seed)
        indices = rng.choice(n_windows, size=n_inject, replace=False)
    else:
        if not 0 <= window_index < n_windows:
            raise ValueError(
                "window_index must be in [0, {}), got {}.".format(n_windows, window_index))
        indices = np.array([window_index])

    for w in indices:
        start = w * features_per_window
        end = start + features_per_window
        block = f1[:, start:end]
        std = block.std(axis=0, ddof=1)
        # Only inject into features with genuine variance (see
        # inject_second_order_joint_leak_with_indices for the rationale); biasing a
        # constant feature would explode Welch's |t| without carrying real leakage.
        median_std = float(np.median(std))
        threshold_std = median_std * 1e-3
        live = std > threshold_std
        if live.any():
            block[:, live] += magnitude * std[live]

    return f1, len(indices)


def inject_second_order_joint_leak_with_indices(f0, f1, window_size, per_feature_sigma=0.05,
                                                fraction=0.1, window_index=None, seed=0):
    """Inject a JOINT-ONLY second-order leak: sub-threshold per feature, visible jointly.

    For a chosen subset of trace windows, adds a small bias to ALL centered-product
    features of the random-group matrix ``f1``. The bias per feature is
    ``per_feature_sigma`` standard deviations of that feature. By construction:

    - Per-feature Welch ``|t|`` stays BELOW the TVLA threshold (``c * sqrt(n/2) < 4.5``
      for ``per_feature_sigma = c`` and ``n`` traces/group) -> the univariate test is
      BLIND to each individual feature.
    - The per-window Hotelling T^2 statistic (Mahalanobis distance
      ``D^2 ~ sum_j (Delta_j / sigma_j)^2 ~ features_per_window * c^2``) exceeds the F
      critical value -> the multivariate test DETECTS the leak.

    This is the canonical "multivariate strictly beats univariate" construction: it is
    the empirical demonstration that the unified multivariate framework gives stronger
    leakage-certification guarantees than per-coordinate univariate TVLA.

    Parameters
    ----------
    f0, f1 : np.ndarray, shape (n0, F) / (n1, F)
        Second-order feature matrices (output of
        :func:`preprocess_second_order_traces`). ``f1`` is modified in place.
    window_size : int
        Trace-point window size; ``F = n_windows * window_size*(window_size+1)//2``.
    per_feature_sigma : float
        Per-feature bias in **standard deviations** of that feature. Must be ``> 0``.
        For balanced ``n`` traces/group, per-feature Welch ``|t| ~ per_feature_sigma *
        sqrt(n/2)``; choose ``per_feature_sigma`` so this stays below 4.5 (e.g.
        ``0.3`` for ``n = 256`` gives ``|t| ~ 3.4``).
    fraction : float
        Fraction of windows to leak (``0 < fraction <= 1``). Ignored when
        ``window_index`` is given.
    window_index : int, optional
        If given, leak exactly this window (must be in ``[0, n_windows)``).
    seed : int
        Seed for the RNG picking which windows are leaked when using ``fraction``.

    Returns
    -------
    (f1, inject_indices)
        ``f1`` (the same array, mutated) and ``inject_indices`` -- a 1-D ``np.ndarray``
        of the window indices that were injected. Use it to run a per-window audit
        (which injected windows did Welch miss but Hotelling catch?).
    """
    features_per_window = window_size * (window_size + 1) // 2
    n_windows = f1.shape[1] // features_per_window
    if per_feature_sigma <= 0:
        raise ValueError("per_feature_sigma must be > 0, got {}.".format(per_feature_sigma))

    if window_index is None:
        if not 0 < fraction <= 1:
            raise ValueError("fraction must be in (0, 1], got {}.".format(fraction))
        n_inject = int(round(n_windows * fraction))
        rng = np.random.RandomState(seed)
        indices = rng.choice(n_windows, size=n_inject, replace=False)
    else:
        if not 0 <= window_index < n_windows:
            raise ValueError(
                "window_index must be in [0, {}), got {}.".format(n_windows, window_index))
        indices = np.array([window_index])

    for w in indices:
        start = w * features_per_window
        end = start + features_per_window
        block = f1[:, start:end]
        std = block.std(axis=0, ddof=1)
        # Only inject into features that have genuine variance. A centered-product
        # feature whose std is ~0 (rank-deficient windows common under ELMO's affine
        # power model) is a *constant* feature -- biasing it would (a) carry no real
        # "leak" information and (b) explode Welch's |t| to ~1e15 (a large mean shift
        # on a zero-variance column -> divide by ~0). The window's median std marks
        # the live scale; features below it are left untouched.
        median_std = float(np.median(std))
        threshold_std = median_std * 1e-3
        live = std > threshold_std
        if live.any():
            block[:, live] += per_feature_sigma * std[live]

    return f1, np.asarray(indices, dtype=np.int64)


# --------------------------------------------------------------------------
# Numerical helpers for the F distribution (no scipy)
# --------------------------------------------------------------------------

def _betacf(a, b, x, itmax=200, eps=3e-12):
    """Continued fraction for the regularized incomplete beta function I_x(a, b).

    Implements ``betacf`` from Numerical Recipes (Press et al.), used by
    `_regularized_beta`. Returns I_x(a, b) * gamma(a+b) / (gamma(a) gamma(b)) * x^a * (1-x)^b
    without the prefactor; see `_regularized_beta`.
    """
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _regularized_beta(a, b, x):
    """Regularized incomplete beta function I_x(a, b) = B(a,b;x) / B(a,b).

    Uses the Numerical Recipes continued fraction, selecting the evaluation
    branch that keeps the continued fraction convergent.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _beta_quantile(a, b, p, itmax=200, tol=1e-14):
    """Quantile of the Beta(a, b) distribution: x with I_x(a, b) = p.

    Bisection on the monotone regularized incomplete beta function.
    """
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1), got {}.".format(p))
    lo, hi = 0.0, 1.0
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        if _regularized_beta(a, b, mid) < p:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)
