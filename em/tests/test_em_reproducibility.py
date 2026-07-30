# ABOUTME: Verify compute_M1_BIC and compute_null_BIC produce identical results
# ABOUTME: when re-run on the same input (deterministic seed from X content).

import numpy as np

from EMfunctions_ASM import compute_M1_BIC, compute_null_BIC


def _make_X(I=40, J=10, p=0.3, seed=123):
    rng = np.random.default_rng(seed)
    X = rng.binomial(1, p, size=(I, J)).astype(float)
    tags = np.array(["noH"] * I)
    return X, tags


def test_compute_M1_BIC_is_reproducible():
    X, tags = _make_X()
    bic_a, theta_a, n_a = compute_M1_BIC(X, tags)
    bic_b, theta_b, n_b = compute_M1_BIC(X, tags)
    assert bic_a == bic_b
    assert n_a == n_b
    np.testing.assert_array_equal(theta_a, theta_b)


def test_compute_null_BIC_is_reproducible():
    X, _ = _make_X()
    a = compute_null_BIC(X)
    b = compute_null_BIC(X)
    assert a == b


def test_M1_BIC_differs_on_different_input():
    X1, tags = _make_X(seed=1)
    X2, _ = _make_X(seed=2)
    _, theta_1, _ = compute_M1_BIC(X1, tags)
    _, theta_2, _ = compute_M1_BIC(X2, tags)
    assert not np.array_equal(theta_1, theta_2)
