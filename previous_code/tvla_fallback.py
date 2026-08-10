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
