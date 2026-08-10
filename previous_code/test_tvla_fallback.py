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
    print("\nAll checks passed!")


if __name__ == "__main__":
    main()
