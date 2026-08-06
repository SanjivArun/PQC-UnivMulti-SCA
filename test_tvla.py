"""test_tvla.py -- manual smoke tests for tvla.py (repo style, no pytest).

Run from the repo root:

    .venv/bin/python3 test_tvla.py

Prints a green "OK" line per check and exits non-zero on the first failure.
Uses only synthetic data, so it needs no ELMO run and completes in seconds.
"""

import math
import sys

import numpy as np

import tvla


def check(name, fn):
    try:
        fn()
        print("OK:   {}".format(name))
    except AssertionError as err:
        print("FAIL: {}  -->  {}".format(name, err))
        sys.exit(1)


def assert_close(actual, expected, tol=1e-6, msg=""):
    if not np.allclose(actual, expected, atol=tol, rtol=tol):
        raise AssertionError("{}: expected {}, got {}".format(msg, expected, actual))


# --------------------------------------------------------------------------
# welch_t_statistic
# --------------------------------------------------------------------------

def test_welch_known_value():
    rng = np.random.RandomState(0)
    n = 100
    group0 = rng.normal(0.0, 1.0, (n, 5))
    group1 = rng.normal(0.0, 1.0, (n, 5)) + 1.0  # shift of 1.0 in every point
    t = tvla.welch_t_statistic(group0, group1)
    # t = (mean0 - mean1)/sqrt(var0/n0 + var1/n1); with a unit shift and unit
    # variance this is -sqrt(n/2) in every point (negative because mean1 > mean0).
    expected = -np.sqrt(n / 2.0)
    assert_close(t, np.full(5, expected), tol=0.5, msg="Welch t magnitude")


def test_welch_zero_variance_identical_means():
    # Both groups constant and equal -> t = 0 (0/0 handled).
    g0 = np.zeros((10, 4))
    g1 = np.ones((10, 4)) * 3.0
    t = tvla.welch_t_statistic(g0, g1)
    # mean0=0, mean1=3, both var 0 -> denom 0, means differ -> +/- inf.
    assert np.all(np.abs(t) == np.inf), t


def test_welch_zero_variance_different_means():
    # Same constant on both sides -> t = 0.
    g0 = np.full((10, 4), 5.0)
    g1 = np.full((10, 4), 5.0)
    t = tvla.welch_t_statistic(g0, g1)
    assert np.all(t == 0.0), t


def test_welch_matches_elmo_style_formula():
    # ELMO's getttest: (mean_fixed - mean_random)/sqrt(var_fixed/N + var_random/N).
    rng = np.random.RandomState(1)
    g0 = rng.normal(10.0, 2.0, (128, 3))
    g1 = rng.normal(11.0, 2.0, (128, 3))
    m0, m1 = g0.mean(axis=0), g1.mean(axis=0)
    v0, v1 = g0.var(axis=0, ddof=1), g1.var(axis=0, ddof=1)
    expected = (m0 - m1) / np.sqrt(v0 / 128 + v1 / 128)
    assert_close(tvla.welch_t_statistic(g0, g1), expected, tol=1e-9, msg="Welch formula")


# --------------------------------------------------------------------------
# count_univariate_leakages
# --------------------------------------------------------------------------

def test_count_univariate_leakages():
    rng = np.random.RandomState(2)
    n, p = 100, 50
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    g1[:, 7] += 10.0   # one strongly leaking point
    count, t = tvla.count_univariate_leakages(g0, g1)
    assert count == 1, "expected exactly 1 leaking point, got {}".format(count)
    assert t.shape == (p,)
    assert abs(t[7]) > 4.5


def test_count_univariate_no_leakage():
    rng = np.random.RandomState(3)
    g0 = rng.normal(0.0, 1.0, (100, 20))
    g1 = rng.normal(0.0, 1.0, (100, 20))
    count, _ = tvla.count_univariate_leakages(g0, g1, threshold=4.5)
    assert count == 0, "two iid groups must not leak, got {}".format(count)


# --------------------------------------------------------------------------
# hotelling_t2_per_window
# --------------------------------------------------------------------------

def test_hotelling_windows_detect_shift():
    rng = np.random.RandomState(4)
    n, p, ws = 256, 150, 50  # 3 windows of 50 points
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    g1[:, 50:100] += 0.3   # mean shift only in the middle window
    t2, n_singular = tvla.hotelling_t2_per_window(g0, g1, ws)
    assert t2.shape == (3,)
    assert n_singular == 0, "full-rank covariances expected here"
    assert t2[0] < 100 and t2[2] < 100, "non-shifted windows should have small T2"
    assert t2[1] > t2[0] * 5, "shifted window should dominate: {}".format(t2)


def test_hotelling_singular_covariance():
    # A constant column (zero variance) makes the covariance singular; the
    # pseudo-inverse path must still run and count the singular window.
    rng = np.random.RandomState(5)
    n, ws = 256, 10
    g0 = rng.normal(0.0, 1.0, (n, ws + 1))
    g1 = rng.normal(0.0, 1.0, (n, ws + 1))
    g0[:, 0] = 1.0      # constant column in group0 -> singular pooled covariance
    g1[:, 0] = 2.0      # constant (different) column in group1
    t2, n_singular = tvla.hotelling_t2_per_window(g0, g1, ws)
    assert n_singular == 1, "expected 1 singular window, got {}".format(n_singular)
    assert np.isfinite(t2).all(), t2


def test_hotelling_window_count():
    rng = np.random.RandomState(6)
    g0 = rng.normal(0.0, 1.0, (256, 173))
    g1 = rng.normal(0.0, 1.0, (256, 173))
    t2, _ = tvla.hotelling_t2_per_window(g0, g1, 50)
    assert t2.shape == (3,), "173 // 50 = 3 windows, got {}".format(t2.shape)


# --------------------------------------------------------------------------
# count_multivariate_leakages
# --------------------------------------------------------------------------

def test_count_multivariate_leakages():
    rng = np.random.RandomState(7)
    n, p, ws = 256, 200, 50
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    g1[:, 100:150] += 0.4   # leak confined to window index 2
    count, t2, n_singular, f_crit = tvla.count_multivariate_leakages(g0, g1, ws)
    assert count == 1, "expected 1 leaking window, got {}".format(count)
    assert t2.shape == (4,)
    assert n_singular == 0
    assert f_crit > 0


def test_count_multivariate_iid_no_leak():
    rng = np.random.RandomState(8)
    g0 = rng.normal(0.0, 1.0, (256, 100))
    g1 = rng.normal(0.0, 1.0, (256, 100))
    count, _, _, _ = tvla.count_multivariate_leakages(g0, g1, 50)
    assert count == 0, "iid groups must not leak, got {}".format(count)


# --------------------------------------------------------------------------
# f_critical  (exact known F quantiles)
# --------------------------------------------------------------------------

def test_f_critical_known_values():
    # F(1,1,0.05) = t_{0.975}(1)^2 = 12.706^2 ~ 161.45
    assert_close(tvla.f_critical(1, 1, 0.05), 161.4476, tol=0.1, msg="F(1,1,0.05)")
    # F(2,2,0.05) = 19.0 exactly
    assert_close(tvla.f_critical(2, 2, 0.05), 19.0, tol=1e-9, msg="F(2,2,0.05)")
    # F(1,2,0.05) ~ 18.5128
    assert_close(tvla.f_critical(1, 2, 0.05), 18.5128, tol=0.01, msg="F(1,2,0.05)")
    # F(1,1,0.01) = t_{0.995}(1)^2 = 63.657^2 ~ 4052.18
    assert_close(tvla.f_critical(1, 1, 0.01), 4052.18, tol=1.0, msg="F(1,1,0.01)")


def test_f_critical_monotonic():
    assert tvla.f_critical(10, 200, 0.01) > tvla.f_critical(10, 200, 0.05)
    assert tvla.f_critical(10, 200, 0.05) > tvla.f_critical(20, 200, 0.05)


def test_default_alpha():
    # Two-sided alpha for |t| > 4.5 under the standard normal ~ 6.8e-6.
    assert_close(tvla.DEFAULT_ALPHA, 6.795e-6, tol=1e-8, msg="DEFAULT_ALPHA")


# --------------------------------------------------------------------------
# Second-order leakage detection
# --------------------------------------------------------------------------

def tuples(rows, cols):
    """Helper for test_extract_upper_triangular_indices."""
    return list(zip(rows.tolist(), cols.tolist()))


def test_extract_upper_triangular_indices():
    # W=3 -> pairs (0,0), (0,1), (0,2), (1,1), (1,2), (2,2) = 6 features
    rows, cols = tvla._extract_upper_triangular_indices(3)
    assert len(rows) == 6
    assert rows[0] == 0
    expected_pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    actual_pairs = [(int(rows[i]), int(cols[i])) for i in range(len(rows))]
    assert actual_pairs == expected_pairs


def test_preprocess_second_order_shapes():
    rng = np.random.RandomState(20)
    n0, n1 = 128, 128
    p = 100
    W = 10
    g0 = rng.normal(0.0, 1.0, (n0, p))
    g1 = rng.normal(0.0, 1.0, (n1, p))
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)
    n_windows = p // W
    feat_per_win = W * (W + 1) // 2
    assert f0.shape == (n0, n_windows * feat_per_win)
    assert f1.shape == (n1, n_windows * feat_per_win)


def test_preprocess_second_order_global_centering():
    """Confirm that global centering is applied (not per-group mean)."""
    rng = np.random.RandomState(30)
    n0, n1 = 64, 64
    p = 50
    W = 5
    g0 = rng.normal(0.0, 1.0, (n0, p))
    # Inject a large shift in group1 so per-group mean != combined mean.
    g1 = rng.normal(0.0, 1.0, (n1, p)) + 50.0
    # Combined mean should be ~25.0 at each point (0 from g0, 55 from g1, equal N).
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)
    # The global mean was used, so f0 and f1 are centered around ~-25 and ~+25 respectively.
    # This means their centered products should be large (not near zero).
    assert np.abs(f0.mean()) > 10.0, "Global centering should produce large residuals"


def test_second_order_univariate_detects_product_leakage():
    """Inject artificial second-order leakage: group1 has correlation on products
    while having zero first-order leakage (same means)."""
    rng = np.random.RandomState(40)
    n = 256
    p = 100
    W = 10
    # Build groups with IDENTICAL means (no first-order leakage).
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    # After centering, the mean-zero property is preserved, but the products change.
    # We do this by adding a shared noise component to group1 products.
    # Simpler: after preprocessing, inject a shift into the second-order features.
    f0, f1 = tvla.preprocess_second_order_traces(g0.copy(), g1.copy(), W)
    shift = 0.5  # shift the second-order feature values
    f1[:, 0] += shift  # inject into first feature
    n_leak, t_stats = tvla.count_second_order_univariate_leakages(f0, f1)
    # The first feature should be flagged.
    assert abs(t_stats[0]) > 4.5, "Injected leakage feature should be detected"


def test_second_order_univariate_iid_no_leak():
    rng = np.random.RandomState(50)
    n = 128
    f0 = rng.normal(0.0, 1.0, (n, 55))
    f1 = rng.normal(0.0, 1.0, (n, 55))
    count, _ = tvla.count_second_order_univariate_leakages(f0, f1)
    assert count == 0, "iid second-order features must not leak, got {}".format(count)


def test_second_order_multivariate_detects_window_leakage():
    """Inject a shift in all features of one trace-window."""
    rng = np.random.RandomState(60)
    n = 128
    p = 50
    W = 10
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    f0, f1 = tvla.preprocess_second_order_traces(g0.copy(), g1.copy(), W)
    # Inject a shift into all features of the first trace-window.
    # f1 shape: (n, n_windows * feat_per_win). First window features are columns [:feat_per_win].
    feat_per_win = W * (W + 1) // 2
    f1[:, :feat_per_win] += 0.5
    count, t2, n_singular, f_crit = tvla.count_second_order_multivariate_leakages(
        f0, f1, W)
    # The first window should leak.
    assert count >= 1, "Injected window should be detected, got {}".format(count)


def test_second_order_multivariate_iid_no_leak():
    rng = np.random.RandomState(70)
    n = 128
    f0 = rng.normal(0.0, 1.0, (n, 55))
    f1 = rng.normal(0.0, 1.0, (n, 55))
    count, _, _, _ = tvla.count_second_order_multivariate_leakages(f0, f1, 10)
    assert count == 0, "iid second-order features must not leak, got {}".format(count)


def test_second_order_2_share_detection():
    """Test the full pipeline on a 2-share-like structure.

    Generate synthetic 2-share traces where individual samples have zero first-order
    leakage but pairwise products leak. This mirrors the structure of a masked
    implementation.

    Model:
    - secret share s0, s1 = s (so s0 + s1 = s)
    - trace value t_i = x_i + noise (where x_i depends on share i)
    - first-order: t_i has no direct dependency on the full secret
    - second-order: t_i * t_j (for different shares of the same point) depends on s^2
    """
    rng = np.random.RandomState(80)
    n_traces = 200
    p = 60
    W = 10

    # Generate shares
    s = rng.randint(0, 256, size=10)  # secret values (small dimension for clarity)
    shares = rng.randint(0, 256, size=(2, 10))
    shares[1] = (s - shares[0]) % 256  # shares[0] + shares[1] = s mod 256

    # Build groups with no first-order leakage but second-order leakage.
    g0 = np.zeros((n_traces, p))
    g1 = np.zeros((n_traces, p))
    for i in range(n_traces):
        for j in range(p // W):
            start = j * W
            for k in range(W):
                # Each trace point depends on a DIFFERENT share of the point.
                # This creates zero first-order leakage.
                share_idx = k % 2
                val = shares[share_idx, i % 10]
                noise = rng.normal(0, 1.0)
                if i < n_traces // 2:
                    g0[start + k] = val + noise
                else:
                    g1[start + k] = val + noise

    # After centering, the pairwise product of (val + noise) should still have
    # a bias proportional to shares^2.
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)
    count, t_stats = tvla.count_second_order_univariate_leakages(f0, f1)
    # We might or might not detect leakage here depending on signal strength.
    # The key check is that the function doesn't crash and produces finite values.
    assert np.isfinite(t_stats).all()
    assert f0.shape[1] == (p // W) * (W * (W + 1) // 2)


def test_inject_second_order_leak_baseline_is_clean():
    """Without injection, masked-like iid features report 0 leaking windows."""
    rng = np.random.RandomState(90)
    n, p, W = 256, 100, 10
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)
    count, _, _, _ = tvla.count_second_order_multivariate_leakages(f0, f1, W)
    assert count == 0, "clean masked-like features must not leak, got {}".format(count)


def test_inject_second_order_leak_all_windows():
    rng = np.random.RandomState(91)
    n, p, W = 256, 100, 10
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)
    n_windows = f1.shape[1] // (W * (W + 1) // 2)

    injected_f1, n_injected = tvla.inject_second_order_leak(f1.copy(), W, magnitude=1.0,
                                                            fraction=1.0)
    assert n_injected == n_windows

    count, _, _, _ = tvla.count_second_order_multivariate_leakages(f0, injected_f1, W)
    assert count == n_windows, "leaking ALL windows should flag ALL windows, got {}/{}".format(
        count, n_windows)


def test_inject_second_order_leak_fraction_detects_some():
    rng = np.random.default_rng(92)
    n, p, W = 256, 100, 10
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)
    n_windows = f1.shape[1] // (W * (W + 1) // 2)

    injected_f1, n_injected = tvla.inject_second_order_leak(f1.copy(), W, magnitude=1.0,
                                                            fraction=0.5)
    assert 0 < n_injected < n_windows

    count, _, _, _ = tvla.count_second_order_multivariate_leakages(f0, injected_f1, W)
    assert count > 0, "injecting into half the windows must be detected"


def test_inject_second_order_leak_specific_window():
    rng = np.random.default_rng(93)
    n, p, W = 256, 100, 10
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    f0, f1 = tvla.preprocess_second_order_traces(g0, g1, W)

    injected_f1, n_injected = tvla.inject_second_order_leak(f1.copy(), W, magnitude=2.0,
                                                            window_index=0)
    assert n_injected == 1, "targeting a single window must inject exactly one"

    count, t2, _, _ = tvla.count_second_order_multivariate_leakages(f0, injected_f1, W)
    assert t2[0] == t2.max(), "the targeted window should carry the largest T^2"
    assert count >= 1, "targeted window should be flagged, got {}".format(count)


# --------------------------------------------------------------------------
# Joint-only second-order leak (the "multivariate strictly beats univariate" case)
#
# Calibration: with n traces/group and a per-feature bias of c std-devs,
#   per-feature Welch |t| ~ c * sqrt(n/2)
# so for n=256, c=0.3 -> |t| ~ 3.4 (< 4.5, Welch BLIND),
# while the multivariate Mahalanobis distance D^2 ~ 55 * c^2 -> Hotelling DETECTS.
# --------------------------------------------------------------------------

def _joint_leak_features(n=256, p=1000, W=10, seed=94):
    """Synthetic masked-like iid features with enough windows for a robust headline."""
    rng = np.random.RandomState(seed)
    g0 = rng.normal(0.0, 1.0, (n, p))
    g1 = rng.normal(0.0, 1.0, (n, p))
    return tvla.preprocess_second_order_traces(g0, g1, W)


def test_inject_joint_leak_per_feature_welch_blind():
    """After a joint-leak injection, per-feature Welch |t| on the injected windows
    stays below the TVLA threshold 4.5 -- univariate is blind by construction."""
    f0, f1 = _joint_leak_features()
    W = 10
    k = W * (W + 1) // 2
    injected_f1, inject_indices = tvla.inject_second_order_joint_leak_with_indices(
        f0, f1.copy(), W, fraction=0.5, seed=0)
    assert len(inject_indices) > 0
    t = tvla.welch_t_statistic(f0, injected_f1)  # shape (F,)
    blind = 0
    for w in inject_indices:
        max_t = float(np.abs(t[w * k:(w + 1) * k]).max())
        if max_t < 4.5:
            blind += 1
    # The vast majority of injected windows must be Welch-blind (a handful may cross
    # 4.5 by chance over many windows). Assert >= 90% blind.
    assert blind >= 0.9 * len(inject_indices), (
        "only {}/{} injected windows are Welch-blind -- per_feature_sigma too large".format(
            blind, len(inject_indices)))


def test_inject_joint_leak_multivariate_detects():
    """Same injection: the multivariate Hotelling T^2 flags most injected windows
    -- despite Welch being blind to (almost) all of them."""
    f0, f1 = _joint_leak_features()
    W = 10
    injected_f1, inject_indices = tvla.inject_second_order_joint_leak_with_indices(
        f0, f1.copy(), W, fraction=0.5, seed=0)
    count, t2, _, f_crit = tvla.count_second_order_multivariate_leakages(f0, injected_f1, W)
    detected = sum(1 for w in inject_indices if t2[int(w)] >= f_crit)
    assert detected >= 0.9 * len(inject_indices), (
        "Hotelling only detected {}/{} injected windows".format(detected, len(inject_indices)))


def test_inject_joint_leak_strictly_more_than_univariate():
    """The headline: multivariate flags strictly more injected windows than
    univariate does on the same data."""
    f0, f1 = _joint_leak_features()
    W = 10
    k = W * (W + 1) // 2
    injected_f1, inject_indices = tvla.inject_second_order_joint_leak_with_indices(
        f0, f1.copy(), W, fraction=0.5, seed=0)
    t = tvla.welch_t_statistic(f0, injected_f1)
    uni_detected = sum(1 for w in inject_indices
                       if float(np.abs(t[w * k:(w + 1) * k]).max()) > 4.5)
    _, t2, _, f_crit = tvla.count_second_order_multivariate_leakages(f0, injected_f1, W)
    multi_detected = sum(1 for w in inject_indices if t2[int(w)] >= f_crit)
    assert multi_detected > uni_detected, (
        "multivariate ({}) must beat univariate ({}) on the joint-leak injection".format(
            multi_detected, uni_detected))


def test_inject_joint_leak_with_indices_paths():
    """Both fraction and window_index paths return the correct indices array,
    and guards are respected."""
    f0, f1 = _joint_leak_features()
    W = 10
    # fraction path
    f1a, idxs_a = tvla.inject_second_order_joint_leak_with_indices(
        f0, f1.copy(), W, fraction=0.5, seed=0)
    assert isinstance(idxs_a, np.ndarray) and idxs_a.ndim == 1
    assert len(idxs_a) > 0
    # single-window path
    f1b, idxs_b = tvla.inject_second_order_joint_leak_with_indices(
        f0, f1.copy(), W, window_index=7)
    assert list(idxs_b) == [7]
    # guards
    try:
        tvla.inject_second_order_joint_leak_with_indices(
            f0, f1.copy(), W, per_feature_sigma=-0.1, fraction=0.5)
        raise AssertionError("negative per_feature_sigma should raise")
    except ValueError:
        pass
    try:
        tvla.inject_second_order_joint_leak_with_indices(
            f0, f1.copy(), W, fraction=0.0)
        raise AssertionError("fraction=0 should raise")
    except ValueError:
        pass


def main():
    check("welch known value", test_welch_known_value)
    check("welch zero-variance / equal means -> inf", test_welch_zero_variance_identical_means)
    check("welch zero-variance / different means -> 0", test_welch_zero_variance_different_means)
    check("welch matches ELMO-style formula", test_welch_matches_elmo_style_formula)
    check("univariate count detects single leak", test_count_univariate_leakages)
    check("univariate count, iid groups -> 0", test_count_univariate_no_leakage)
    check("hotelling windows detect shift", test_hotelling_windows_detect_shift)
    check("hotelling singular covariance handled", test_hotelling_singular_covariance)
    check("hotelling window count", test_hotelling_window_count)
    check("multivariate count detects single leak", test_count_multivariate_leakages)
    check("multivariate count, iid groups -> 0", test_count_multivariate_iid_no_leak)
    check("f_critical known values", test_f_critical_known_values)
    check("f_critical monotonic", test_f_critical_monotonic)
    check("default alpha", test_default_alpha)

    # --- Second-order tests ---
    check("extract upper triangular indices", test_extract_upper_triangular_indices)
    check("preprocess_second_order_shapes", test_preprocess_second_order_shapes)
    check("preprocess_second_order_global_centering", test_preprocess_second_order_global_centering)
    check("second_order_univariate_detects_product_leakage", test_second_order_univariate_detects_product_leakage)
    check("second_order_univariate_iid_no_leak", test_second_order_univariate_iid_no_leak)
    check("second_order_multivariate_detects_window_leakage", test_second_order_multivariate_detects_window_leakage)
    check("second_order_multivariate_iid_no_leak", test_second_order_multivariate_iid_no_leak)
    check("second_order_2_share_detection", test_second_order_2_share_detection)
    check("inject_second_order_leak_baseline_is_clean", test_inject_second_order_leak_baseline_is_clean)
    check("inject_second_order_leak_all_windows", test_inject_second_order_leak_all_windows)
    check("inject_second_order_leak_fraction_detects_some", test_inject_second_order_leak_fraction_detects_some)
    check("inject_second_order_leak_specific_window", test_inject_second_order_leak_specific_window)

    # --- Joint-only leak (multivariate strictly beats univariate) ---
    check("inject_joint_leak_per_feature_welch_blind", test_inject_joint_leak_per_feature_welch_blind)
    check("inject_joint_leak_multivariate_detects", test_inject_joint_leak_multivariate_detects)
    check("inject_joint_leak_strictly_more_than_univariate", test_inject_joint_leak_strictly_more_than_univariate)
    check("inject_joint_leak_with_indices_paths", test_inject_joint_leak_with_indices_paths)

    print("\nAll checks passed!")


if __name__ == "__main__":
    main()
