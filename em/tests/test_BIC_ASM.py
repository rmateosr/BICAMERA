# ABOUTME: Regression test — Shared (M0) and ASM (M1) BIC parameter counts (J and 2J).
# ABOUTME: Exercises the BIC formulas through compute_null_BIC and compute_M1_BIC.
import numpy as np

import EMfunctions_ASM as em


def test_null_bic_parameter_count_is_J():
    """Shared model: theta_j shared across all reads. p = J."""
    rng = np.random.RandomState(0)
    J = 5
    X = (rng.rand(30, J) > 0.5).astype(float)

    bic = em.compute_null_BIC(X)

    # Recompute the expected BIC by hand.
    total_meth = X.sum(axis=0)
    total_obs = float(X.shape[0])
    theta = total_meth / (total_obs + em.PSEUDO)
    log_theta = np.log(theta + em.PSEUDO)
    log_1mtheta = np.log(1 - theta + em.PSEUDO)
    ll = float(np.sum(X * log_theta) + np.sum((1 - X) * log_1mtheta))
    expected = J * np.log(X.shape[0]) - 2 * ll

    np.testing.assert_allclose(bic, expected, atol=1e-6)


def test_M1_bic_parameter_count_is_2J():
    """ASM model under the diploid constraint: per-allele theta at k=0, p = 2J."""
    rng = np.random.RandomState(1)
    J = 5
    X = (rng.rand(30, J) > 0.5).astype(float)
    tags = np.array(["H1"] * 10 + ["H2"] * 10 + ["noH"] * 10)

    bic, theta, _ = em.compute_M1_BIC(X, tags)

    # Recompute the expected BIC from the EM's terminal (theta, pi).
    pi = np.zeros((2, 2))
    pi[:, 0] = 0.5
    ll = em._observed_loglikelihood(X, tags, theta, pi)
    expected = 2 * J * np.log(X.shape[0]) - 2 * ll

    np.testing.assert_allclose(bic, expected, atol=1e-6)


def test_bic_formula_strings_present_in_source():
    """Cheap regression check: parameter counts (J and 2J) appear in source."""
    import inspect

    src_null = inspect.getsource(em.compute_null_BIC)
    src_M1 = inspect.getsource(em.compute_M1_BIC)
    assert "p_null = J" in src_null
    assert "p_M1 = 2 * J" in src_M1
