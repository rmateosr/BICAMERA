# ABOUTME: Unit tests for doM_theta on the ASM model.
# ABOUTME: Confirms weighted average semantics and per-(h, k) independence.
import numpy as np

from EMfunctions_ASM import doM_theta, PSEUDO


def test_weighted_average_per_h_k():
    """doM_theta returns the gamma-weighted average of X per (j, h, k)."""
    # 4 reads, 3 CpGs. Hand-built so each (h, k) cell is dominated by one read.
    X = np.array([
        [1.0, 1.0, 0.0],   # read 0 — fully methylated at j=0,1
        [0.0, 0.0, 1.0],   # read 1
        [1.0, 0.0, 1.0],   # read 2
        [0.0, 1.0, 0.0],   # read 3
    ])
    gamma = np.zeros((4, 2, 2))
    gamma[0, 0, 0] = 1.0  # read 0 -> (h=0, k=0)
    gamma[1, 1, 0] = 1.0  # read 1 -> (h=1, k=0)
    gamma[2, 0, 1] = 1.0  # read 2 -> (h=0, k=1)
    gamma[3, 1, 1] = 1.0  # read 3 -> (h=1, k=1)

    theta = doM_theta(X, gamma)

    np.testing.assert_allclose(theta[:, 0, 0], X[0])
    np.testing.assert_allclose(theta[:, 1, 0], X[1])
    np.testing.assert_allclose(theta[:, 0, 1], X[2])
    np.testing.assert_allclose(theta[:, 1, 1], X[3])


def test_nan_sites_excluded_from_update():
    """NaN entries in X should contribute to neither numerator nor denominator.

    With one read whose values at site j=1 are NaN, the (h=0, k=0) theta at j=1
    should fall back to the prior default 0.5 (denominator = 0).
    """
    X = np.array([
        [1.0, np.nan, 0.0],
    ])
    gamma = np.zeros((1, 2, 2))
    gamma[0, 0, 0] = 1.0

    theta = doM_theta(X, gamma)

    np.testing.assert_allclose(theta[0, 0, 0], 1.0)
    np.testing.assert_allclose(theta[1, 0, 0], 0.5)  # NaN site -> default
    np.testing.assert_allclose(theta[2, 0, 0], 0.0)


def test_zero_gamma_keeps_default():
    """A (h, k) cell with zero posterior mass everywhere keeps the default 0.5."""
    X = np.ones((3, 4))
    gamma = np.zeros((3, 2, 2))
    gamma[:, 0, 0] = 1.0  # all mass on (h=0, k=0)

    theta = doM_theta(X, gamma)

    np.testing.assert_allclose(theta[:, 0, 0], 1.0)
    # The other three cells should stay at the 0.5 default — no posterior mass
    # ever placed there, denominator <= PSEUDO.
    np.testing.assert_allclose(theta[:, 1, 0], 0.5)
    np.testing.assert_allclose(theta[:, 0, 1], 0.5)
    np.testing.assert_allclose(theta[:, 1, 1], 0.5)
