# ABOUTME: EM + BIC for the static Shared (M0) and ASM (M1) models on a single sample.
# ABOUTME: ASM uses the diploid per-allele marginal constraint pi[h, 0] = 0.5 (k=1 disallowed).

import hashlib
import os
import re

import numpy as np
import pandas as pd


def _seed_from_X(X):
    """Deterministic 32-bit seed derived from the read matrix.

    Same input matrix -> same seed -> same random theta init -> same EM
    trajectory. Used by compute_M1_BIC to make per-window BIC values
    reproducible across runs.
    """
    h = hashlib.md5()
    h.update(np.ascontiguousarray(X).tobytes())
    return int(h.hexdigest()[:8], 16)


np.seterr(all='raise')

# Debug trace sink for the EM-tol calibration (off in production). When set to a
# list, compute_M1_BIC appends each iteration's observed log-likelihood to it;
# left None during pipeline runs so the convergence loop is byte-identical and
# pays only one `is not None` check per iteration. Set/reset by the calibration
# harness only — never by the driver or SGE scripts.
_M1_LL_TRACE = None

PSEUDO = 1e-10
# Minimum reads per window. Below this the EM is degenerate: with N=1 the BIC
# penalty p*log(N) is zero, so M0 and M1 collapse to the same score and ~1e-9
# floating-point noise in M1's LL decides the winner — yielding spurious ASM
# calls. DEPTH=2 makes the BIC parsimony bonus strictly positive (J*log(2) > 0)
# so M0 wins ties on shape. Higher thresholds (e.g. 5, 10) are a sensitivity
# tuning choice for downstream callers; the EM math is well-defined for N>=2.
DEPTH = 2

# Default EM restart count for the ASM (M1) BIC. Validated on the 2026-05-18
# fine restart sweep against the time-emergent code path: r=5 is where the
# restart-driven BIC drop plateaus. Callers can pass n_restarts=1 to reproduce
# the pre-restart baseline; r=0 reuses the deterministic _seed_from_X seed.
DEFAULT_RESTARTS = 5

# Operational default for the sliding-window step (CpGs advanced between
# consecutive evaluated windows), used by the CLI driver / pipeline. The library
# function below defaults to 1 (evaluate every in-region CpG) so unit tests that
# assert the exact window-start set stay valid; the pipeline overrides to 5 for a
# ~5x faster, coarser scan.
DEFAULT_WINDOW_STEP = 5

# Default absolute convergence tolerance for the ASM (M1) EM: the loop stops when
# abs(ll - previous_ll) < tol. Loosened from 1e-8 to 1e-4 after calibration showed
# 1e-4 to be the loosest tolerance that changes no result: 0 of 73,128 BICresult
# flips on a whole-chromosome HiFi baseline, and bit-identical sensitivity / FPR /
# F1 across the synthetic truth grid including low depth — while cutting EM
# iterations and wall-clock to ~45% of 1e-8. Single fixed absolute value by design.
DEFAULT_EM_TOL = 1e-4


def _m1_best_of_n(X_window, tags_window, n_restarts, tol=DEFAULT_EM_TOL):
    """Best-of-n compute_M1_BIC; monkeypatches _seed_from_X."""
    original = _seed_from_X
    g = globals()
    best = np.inf
    best_iters = 0
    best_theta = None
    try:
        for r in range(n_restarts):
            g["_seed_from_X"] = (
                lambda X, _r=r: (original(X) + _r * 1_000_003) & 0xFFFFFFFF
            )
            bic, theta, n_iters = compute_M1_BIC(X_window, tags_window, tol=tol)
            if bic < best:
                best = bic
                best_iters = n_iters
                best_theta = theta
    finally:
        g["_seed_from_X"] = original
    return float(best), best_theta, best_iters


def Reshape2Matrix(binfullset, coords_array):
    """Convert read-level data to a binary matrix (n_reads x J).

    Each row is a read; each column is a CpG site in coords_array.
    Entries are 0 or 1 (methylation status), NaN if the read doesn't cover that site.

    Parameters
    ----------
    binfullset : DataFrame
        Must have columns 'startcoord' (comma-separated positions) and 'status' (comma-separated 0/1).
    coords_array : ndarray of int
        Sorted CpG positions defining the columns.

    Returns
    -------
    matrix : ndarray, shape (n_reads, len(coords_array))
    """
    nrows = binfullset.shape[0]
    ncols = len(coords_array)
    matrix = np.full((nrows, ncols), np.nan)

    for i in range(nrows):
        row = binfullset.iloc[i]
        read_coords = np.array([int(float(c)) for c in str(row['startcoord']).split(',')])
        read_values = np.array([float(v) for v in str(row['status']).split(',')])

        _, unique_idx = np.unique(read_coords, return_index=True)
        read_coords = read_coords[unique_idx]
        read_values = read_values[unique_idx]

        coord_in_set = np.isin(coords_array, read_coords)
        read_in_set = np.isin(read_coords, coords_array)

        matrix[i, np.where(coord_in_set)[0]] = read_values[read_in_set]

    return matrix


_REGION_FILENAME_RE = re.compile(r"methylationfraction_(\d+)_(\d+)_\.tsv$")


def parse_region_bounds(path):
    """Parse the filename of a per-region TSV into (regionstart, regionend).

    Expected name shape: methylationfraction_<start>_<end>_.tsv. These bounds
    are the strict CpG-coord window the EM driver evaluates inside this
    region (see EMBIC_bin_path_ASM's claim-mask). Hanging CpGs from reads
    that fall outside [start, end] stay in the file (to keep windows at the
    boundary fully supported) but are NOT emitted as coord_evaluated rows.
    """
    fname = os.path.basename(path)
    m = _REGION_FILENAME_RE.match(fname)
    if not m:
        raise ValueError(
            f"cannot parse region bounds from filename: {fname!r} "
            "(expected methylationfraction_<start>_<end>_.tsv)"
        )
    return int(m.group(1)), int(m.group(2))


_BIC_SUMMARY_COLUMNS = (
    "coord_evaluated",
    "window_last_coord",
    "BICsinglecomp",
    "BICmiddlecomp",
    "BICresult",
    "BIC_winner",
    "windows_size",
    "total_reads",
    "em_iterations",
    "mean_theta_diff",
    "window_cpgs",
)


def _empty_bic_summary():
    """Empty BIC summary with the canonical schema."""
    return pd.DataFrame({col: [] for col in _BIC_SUMMARY_COLUMNS})


def EMBIC_bin_path_ASM(path, windowsize=10, asm_restarts=DEFAULT_RESTARTS, window_step=1,
                       tol=DEFAULT_EM_TOL):
    """Run the ASM pipeline on a single region TSV.

    Each region is standalone and non-overlapping (the methylation-guided
    potential-ASM-region route). If the region has >= `windowsize` in-region CpGs,
    every in-region CpG that can start a full W-window is evaluated as a
    window-start (kept every `window_step`-th) — the region's right edge is not
    dropped. If the region has FEWER than `windowsize` in-region CpGs, a single
    window spanning the whole region is evaluated (adaptive window) instead of
    skipping it; sub-3-CpG regions are pruned upstream at the bucket step.
    A region file may still carry hanging CpGs from reads that spill past
    [regionstart, regionend]; those are used for window support but are not
    emitted as coord_evaluated rows.

    Right-edge tail rescue (always on). When `window_step` > 1, the strided
    window-starts leave the last (n_in - windowsize) % window_step in-region CpGs
    of the region in no window — they fall off the right edge and, because regions
    are standalone (the next region starts a fresh slide), are never analyzed. A
    final full-W window is therefore anchored to END on the region's last CpG so
    that orphan tail is always scored in a canonical windowsize-CpG window. Every
    window stays exactly `windowsize` CpGs wide, so the BIC parsimony differential
    J*log(n) is constant across windows and per-window ASM decisions stay
    comparable (the reason a constant-W anchored window is used rather than a
    shorter variable-size tail window).

    Parameters
    ----------
    path : str
        Per-region TSV path produced by bucket_reads_into_regions.py.
        Filename MUST match methylationfraction_<start>_<end>_.tsv.
    windowsize : int
        Number of CpG sites per sliding window.
    asm_restarts : int
        EM restarts for the ASM (M1) BIC; keep min BIC across restarts.
        Default 5; pass 1 to reproduce the pre-restart baseline.
    window_step : int
        CpGs advanced between consecutive evaluated window-starts. Default 1
        (evaluate every in-region CpG). The pipeline passes 5 for a coarser,
        faster scan (see DEFAULT_WINDOW_STEP). Must be >= 1.
    tol : float
        Absolute convergence tolerance for the ASM (M1) EM (stops when
        abs(ll - previous_ll) < tol). Default DEFAULT_EM_TOL (1e-4).

    Returns
    -------
    BIC_summary : DataFrame
    """
    region_start, region_end = parse_region_bounds(path)

    df = pd.read_csv(path, sep='\t')

    coords_str = ','.join(df['startcoord'].astype(str))
    coords_array = np.array(sorted({int(float(c)) for c in coords_str.split(',')}))

    X_full = Reshape2Matrix(df, coords_array)
    tags_full = np.array(df['haplotype'])

    valid = ~np.isnan(X_full).all(axis=1)
    X_full = X_full[valid]
    tags_full = tags_full[valid]

    J_full = len(coords_array)

    # Indices of CpGs that fall inside the filename-declared region.
    in_filename = np.where(
        (coords_array >= region_start) & (coords_array <= region_end)
    )[0]
    if window_step < 1:
        raise ValueError(f"window_step must be >= 1, got {window_step}")
    n_in = in_filename.size
    if n_in == 0:
        return _empty_bic_summary()
    # Adaptive window. Regions are standalone and non-overlapping. If the region
    # has fewer than `windowsize` in-region CpGs, evaluate ONE window spanning the
    # whole region (eff_W = n_in) rather than skipping it — small regions are no
    # longer silently dropped (sub-3-CpG regions are pruned upstream at the bucket
    # step). Otherwise slide a full `windowsize`-wide window, dropping starts whose
    # window would run past coords_array and keeping every window_step-th start.
    # In-region CpG indices are contiguous, so [::window_step] advances exactly
    # window_step CpGs between evaluated windows.
    if n_in < windowsize:
        eff_W = int(n_in)
        claim_indices = in_filename[:1]
    else:
        eff_W = int(windowsize)
        claim_indices = in_filename[in_filename + windowsize <= J_full][::window_step]
        # Right-edge tail rescue: append a final full-W window anchored to END on
        # the last in-region CpG (start = in_filename[-1] - windowsize + 1) so the
        # (n_in - windowsize) % window_step CpGs that striding leaves in no window
        # are still scored. Without this, at window_step=5 ~78% of >=W-CpG regions
        # orphan 1-4 trailing CpGs, dropping real signal at every region right edge
        # (on one HiFi chromosome this recovered 18,634 CpGs / 7,604 windows / 712
        # ASM calls). Added only when the last strided window does not
        # already reach the last in-region CpG. anchor stays in-region (start >=
        # in_filename[0] since n_in >= windowsize) and ends on in_filename[-1]
        # (<= region_end), so no coord_evaluated leaks past the edge.
        if claim_indices.size:
            anchor = in_filename[-1] - windowsize + 1
            if anchor > claim_indices[-1]:
                claim_indices = np.append(claim_indices, anchor)

    n_windows = claim_indices.size
    if n_windows == 0:
        return _empty_bic_summary()

    BICsinglecomp = np.full(n_windows, np.nan)
    BICmiddlecomp = np.full(n_windows, np.nan)
    em_iters = np.full(n_windows, 0, dtype=int)
    theta_diff = np.full(n_windows, np.nan)
    total_reads_per_window = np.zeros(n_windows, dtype=int)

    for out_idx, cont in enumerate(claim_indices):
        X_w = X_full[:, cont:cont + eff_W]
        tags_w = tags_full

        keep = ~np.isnan(X_w).all(axis=1)
        X_window = X_w[keep]
        tags_window = tags_w[keep]
        n_reads = X_window.shape[0]
        total_reads_per_window[out_idx] = n_reads

        if n_reads < DEPTH:
            # Sentinels: argmin of (0, 1) = 0 -> winner = Shared.
            BICsinglecomp[out_idx] = 0
            BICmiddlecomp[out_idx] = 1
            continue

        bic_null, bic_M1, n_iters, mtd = _compute_window_BIC(
            X_window, tags_window, asm_restarts, tol=tol,
        )
        BICsinglecomp[out_idx] = bic_null
        BICmiddlecomp[out_idx] = bic_M1
        em_iters[out_idx] = n_iters
        theta_diff[out_idx] = mtd

    BIC_winner = np.argmin(
        np.stack([BICsinglecomp, BICmiddlecomp], axis=1), axis=1
    )

    # coord_evaluated is the window's FIRST CpG (its start). window_last_coord is
    # the window's last CpG: coords_array[start + eff_W - 1]. For full windows this
    # can sit past region_end (the window extends into hanging CpGs kept for support);
    # for adaptive sub-W regions the window spans the whole region. window_cpgs is the
    # comma-joined list of all eff_W CpG coords in the window (length == windows_size).
    window_last_coord = coords_array[claim_indices + eff_W - 1]
    window_cpgs = [
        ",".join(str(int(c)) for c in coords_array[ci:ci + eff_W])
        for ci in claim_indices
    ]

    BIC_summary = pd.DataFrame({
        "coord_evaluated": coords_array[claim_indices],
        "window_last_coord": window_last_coord,
        "BICsinglecomp": BICsinglecomp,
        "BICmiddlecomp": BICmiddlecomp,
        # Semantics: 1 if ASM (M1) beats Shared (M0). Recovers the original
        # cross-sectional BIC ASM convention (M0 vs M1).
        "BICresult": (BICmiddlecomp < BICsinglecomp).astype(int),
        "BIC_winner": BIC_winner,
        "windows_size": eff_W,
        "total_reads": total_reads_per_window,
        "em_iterations": em_iters,
        "mean_theta_diff": theta_diff,
        "window_cpgs": window_cpgs,
    })
    return BIC_summary


def _compute_window_BIC(X_window, tags_window, asm_restarts=DEFAULT_RESTARTS,
                        tol=DEFAULT_EM_TOL):
    """Run EM and compute BIC for a single window.

    Returns (bic_null, bic_M1, n_iters, mean_theta_diff).

    The ASM (M1) EM is run `asm_restarts` times, keeping the minimum BIC
    across restarts. n_iters and mean_theta_diff come from the restart with
    min BIC.
    """
    bic_null = compute_null_BIC(X_window)
    bic_M1, theta, n_iters = _m1_best_of_n(X_window, tags_window, asm_restarts, tol=tol)

    # mean_theta_diff: average |theta_h0 - theta_h1| across CpGs. Measures
    # how different the per-allele methylation patterns are at k=0 — the
    # primary ASM signal under the diploid constraint.
    # theta shape: (J, 2, 2) — theta[j, h, k]. Only k=0 is used here; k=1
    # is unconstrained noise (pi[:,1] is pinned at 0).
    mean_theta_diff = float(np.mean(np.abs(theta[:, 1, 0] - theta[:, 0, 0])))

    return float(bic_null), float(bic_M1), n_iters, mean_theta_diff


def doE(X, theta, pi, tags):
    """E-step.

    Parameters
    ----------
    X : ndarray, shape (I, J)
        Binary methylation matrix. NaN for missing.
    theta : ndarray, shape (J, 2, 2)
        Current theta estimates. theta[j, h, k] = P(methylated at site j | allele h, model k).
    pi : ndarray, shape (2, 2)
        Current pi estimates. pi[h, k] = P(allele h, model k).
    tags : ndarray of str, shape (I,)
        Haplotype tags: "H1", "H2", or "noH".

    Returns
    -------
    gamma : ndarray, shape (I, 2, 2)
        Posterior probabilities. gamma[i, h, k] = P(z_{i,h,k}=1 | X, theta, pi).
    """
    I, J = X.shape
    gamma = np.zeros((I, 2, 2))

    # NaN -> contribute 0 to log-likelihood via paired methylated/unmethylated masks.
    xformeth = X.copy()
    xformeth[np.isnan(xformeth)] = 0

    xforunmeth = X.copy()
    xforunmeth[np.isnan(xforunmeth)] = 1
    xforunmeth = 1 - xforunmeth

    log_gamma_prime = np.full((I, 2, 2), -np.inf)

    for h in range(2):
        for k in range(2):
            if pi[h, k] < PSEUDO:
                continue

            log_theta = np.log(theta[:, h, k] + PSEUDO)
            log_1_minus_theta = np.log(1 - theta[:, h, k] + PSEUDO)

            ll_per_read = xformeth @ log_theta + xforunmeth @ log_1_minus_theta

            log_gamma_prime[:, h, k] = np.log(pi[h, k] + PSEUDO) + ll_per_read

    h1_mask = tags == "H1"
    h2_mask = tags == "H2"
    noh_mask = ~(h1_mask | h2_mask)

    if np.any(h1_mask):
        gamma[h1_mask, 0, :] = _normalize_log_batch(
            log_gamma_prime[h1_mask, 0, :], axis=-1
        )

    if np.any(h2_mask):
        gamma[h2_mask, 1, :] = _normalize_log_batch(
            log_gamma_prime[h2_mask, 1, :], axis=-1
        )

    if np.any(noh_mask):
        n_noh = int(noh_mask.sum())
        flat = log_gamma_prime[noh_mask].reshape(n_noh, 4)
        gamma[noh_mask] = _normalize_log_batch(flat, axis=-1).reshape(n_noh, 2, 2)

    return gamma


def _normalize_log_batch(log_values, axis=-1):
    """Vectorized softmax over `axis` with -inf masking.

    Slices that are all -inf (or whose exp-shifted sum is 0) return uniform;
    otherwise returns normalized probabilities summing to 1 along `axis`.
    """
    finite_mask = np.isfinite(log_values)
    any_finite = np.any(finite_mask, axis=axis, keepdims=True)

    safe_for_max = np.where(finite_mask, log_values, -np.inf)
    L_star = np.max(safe_for_max, axis=axis, keepdims=True)
    L_star = np.where(any_finite, L_star, 0.0)

    shifted = np.clip(log_values - L_star, -700, None)
    exp_shifted = np.where(finite_mask, np.exp(shifted), 0.0)
    total = exp_shifted.sum(axis=axis, keepdims=True)

    n_along = log_values.shape[axis]
    uniform = np.full_like(log_values, 1.0 / n_along)
    safe_total = np.where(total > 0, total, 1.0)
    normalized = exp_shifted / safe_total
    degenerate = (total == 0) | (~any_finite)
    return np.where(degenerate, uniform, normalized)


def _logsumexp_batch(log_values, axis=-1):
    """Vectorized log-sum-exp over `axis` with -inf masking."""
    finite_mask = np.isfinite(log_values)
    any_finite = np.any(finite_mask, axis=axis, keepdims=True)

    safe_for_max = np.where(finite_mask, log_values, -np.inf)
    L_star = np.max(safe_for_max, axis=axis, keepdims=True)
    L_star = np.where(any_finite, L_star, 0.0)

    shifted = np.clip(log_values - L_star, -700, None)
    exp_shifted = np.where(finite_mask, np.exp(shifted), 0.0)
    total = exp_shifted.sum(axis=axis, keepdims=True)

    result = np.squeeze(L_star + np.log(total + PSEUDO), axis=axis)
    any_finite_squeezed = np.squeeze(any_finite, axis=axis)
    return np.where(any_finite_squeezed, result, 0.0)


def doM_theta(X, gamma):
    """M-step for theta.

    theta[j, h, k] = sum_i gamma[i,h,k] * x_{i,j} / sum_i gamma[i,h,k] * observed_{i,j}

    Parameters
    ----------
    X : ndarray, shape (I, J)
    gamma : ndarray, shape (I, 2, 2)

    Returns
    -------
    theta : ndarray, shape (J, 2, 2)
    """
    J = X.shape[1]
    theta = np.full((J, 2, 2), 0.5)

    observed = ~np.isnan(X)
    X_clean = X.copy()
    X_clean[np.isnan(X_clean)] = 0
    observed_f = observed.astype(float)

    for h in range(2):
        for k in range(2):
            g = gamma[:, h, k]
            numerator = (g[:, np.newaxis] * X_clean).sum(axis=0)
            denominator = (g[:, np.newaxis] * observed_f).sum(axis=0)
            has_data = denominator > PSEUDO
            theta[has_data, h, k] = numerator[has_data] / denominator[has_data]

    return theta


def _observed_loglikelihood(X, tags, theta, pi):
    """Observed-data log-likelihood for the diploid-aware mixture.

    For tagged reads (A1, A2) we *condition* on the haplotype tag and drop
    the log(pi) prior. For untagged reads (A0) the latent allele is
    marginalised out, so log(pi) stays inside the logsumexp. This keeps
    the LL symmetric with M0's marginal LL of X for the BIC comparison.

    Used by compute_M1_BIC to monitor convergence and to score the BIC.
    """
    I, J = X.shape

    xformeth = X.copy()
    xformeth[np.isnan(xformeth)] = 0
    xforunmeth = X.copy()
    xforunmeth[np.isnan(xforunmeth)] = 1
    xforunmeth = 1 - xforunmeth

    # ll_per_read_table[i, h, k] = data-only log-likelihood of read i under
    # (h, k), with no log(pi) prefactor. Slots where pi[h, k] is zero stay
    # at -inf so they are masked out of any later logsumexp.
    ll_per_read_table = np.full((I, 2, 2), -np.inf)
    log_pi_terms = np.full((2, 2), -np.inf)
    for h in range(2):
        for k in range(2):
            if pi[h, k] < PSEUDO:
                continue
            log_theta = np.log(theta[:, h, k] + PSEUDO)
            log_1_minus_theta = np.log(1 - theta[:, h, k] + PSEUDO)
            ll_per_read_table[:, h, k] = xformeth @ log_theta + xforunmeth @ log_1_minus_theta
            log_pi_terms[h, k] = np.log(pi[h, k] + PSEUDO)

    h1_mask = tags == "H1"
    h2_mask = tags == "H2"
    noh_mask = ~(h1_mask | h2_mask)

    ll = 0.0
    # A1: h=0 is observed; only k=0 is active under the diploid constraint.
    if np.any(h1_mask):
        ll += float(ll_per_read_table[h1_mask, 0, 0].sum())
    # A2: h=1 is observed; only k=0 is active.
    if np.any(h2_mask):
        ll += float(ll_per_read_table[h2_mask, 1, 0].sum())
    # A0: marginalise over (h, k) using the pi prior.
    if np.any(noh_mask):
        n_noh = int(noh_mask.sum())
        log_terms_noh = ll_per_read_table[noh_mask] + log_pi_terms[None, :, :]
        flat = log_terms_noh.reshape(n_noh, 4)
        ll += float(_logsumexp_batch(flat, axis=-1).sum())

    return ll


def compute_null_BIC(X):
    """BIC for the Shared model: a single theta_j shared across all reads.

    Parameters
    ----------
    X : ndarray, shape (I, J)

    Returns
    -------
    BIC_null : float
    """
    J = X.shape[1]

    observed = ~np.isnan(X)
    X_clean = X.copy()
    X_clean[np.isnan(X_clean)] = 0
    total_meth = X_clean.sum(axis=0)
    total_obs = observed.astype(float).sum(axis=0)
    total_reads = X.shape[0]

    theta_null = total_meth / (total_obs + PSEUDO)

    X_unmeth = observed.astype(float) - X_clean
    log_theta = np.log(theta_null + PSEUDO)
    log_1_minus_theta = np.log(1 - theta_null + PSEUDO)

    ll_null = float(np.sum(X_clean * log_theta) + np.sum(X_unmeth * log_1_minus_theta))

    p_null = J
    BIC_null = p_null * np.log(total_reads) - 2 * ll_null
    return float(BIC_null)


def compute_M1_BIC(X, tags, maxIter=1000, tol=DEFAULT_EM_TOL):
    """BIC for the ASM model under the diploid per-allele marginal constraint.

    pi[h, 0] = 0.5, pi[h, 1] = 0 pinned for both alleles, so the altered
    component is disallowed and components collapse to a single per-allele
    state at k=0. theta[:, :, 1] is unused. Parameter count is 2J (per-allele
    theta at k=0).

    Parameters
    ----------
    X : ndarray, shape (I, J)
    tags : ndarray of str, shape (I,)
    maxIter : int
    tol : float

    Returns
    -------
    BIC_M1 : float
    theta : ndarray, shape (J, 2, 2)
        Final theta — caller uses theta[:, :, 0] to compute mean_theta_diff.
    n_iters : int
    """
    J = X.shape[1]

    np.random.seed(_seed_from_X(X))
    theta = np.random.uniform(0.01, 0.99, size=(J, 2, 2))

    pi = np.zeros((2, 2))
    pi[:, 0] = 0.5

    previous_ll = -np.inf
    ll = -np.inf
    iteration = 0
    for iteration in range(maxIter):
        gamma = doE(X, theta, pi, tags)
        theta = doM_theta(X, gamma)
        ll = _observed_loglikelihood(X, tags, theta, pi)
        if _M1_LL_TRACE is not None:
            _M1_LL_TRACE.append(ll)
        if abs(ll - previous_ll) < tol:
            break
        previous_ll = ll

    n_total = X.shape[0]
    p_M1 = 2 * J
    BIC_M1 = p_M1 * np.log(n_total) - 2 * ll
    return float(BIC_M1), theta, iteration + 1
