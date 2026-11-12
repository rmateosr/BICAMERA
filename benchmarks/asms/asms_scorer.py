# ABOUTME: Faithful Python reimplementation of ASMS's cluster+filter decision on a read x CpG matrix.
# ABOUTME: K=2 multivariate-Bernoulli EM (10 iters, free pi, single init) + heuristics + 5000-shuffle permutation gate.
#
# Adapted from ASMS by E. Raineri, M. A. Esteve Marco, E. Casals and A. Esteve Codina:
# "asms: finding allele-specific methylation in human genomes without phasing",
# NAR Genomics and Bioinformatics 8(3):lqag072 (2026), doi:10.1093/nargab/lqag072.
# Source: ecmra/asms Zenodo release, doi:10.5281/zenodo.20323001 (commit c7178db),
# used under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
# Changes: the `cluster` and `filter` steps are rewritten in Python; `scan` and BAM
# reading are not ported; the permutation shuffles come from numpy, so p-values
# match the Rust tool only up to Monte-Carlo error.
"""ASMS (Raineri et al. 2026, NAR Genom Bioinform, doi:10.1093/nargab/lqag072) as a sweep scorer.

Reimplements the *decision logic* of `asms cluster` + `asms filter` so it runs
on the same (X) read x CpG matrices the BIC sweep already generates — a true
apples-to-apples comparison on identical reads. Verified against the ecmra/asms
`release` source (commit c7178db) and validated against the real binary via the
matrix-file + `asms pval` path (see validate_asms_scorer.py).

Faithful to the source (NOT the paper prose) on every detail that affects output:

  * EM: K=2 multivariate Bernoulli, free pi re-estimated each iter, exactly
    `maxit` iterations (default 10), NO convergence check, single random init
    (one restart). (cluster.rs)
  * M-step denominator is the GLOBAL cluster responsibility pop[k] = sum_i
    gamma[i,k], not the per-position observed responsibility. (cluster.rs update_mu)
  * log-Bernoulli term has NO pseudocount; special-cases x==mu==0/1 -> 0. (prob.rs)
  * Missing CpGs (NaN here, -1 in Rust) are skipped in the read likelihood and
    in the cluster-mean numerator. (prob.rs logmvbern, cluster.rs update_mu)
  * absmdiff = MEDIAN across CpGs of per-CpG |m1[j]-m0[j]|, dropping CpGs with no
    reads in a cluster (NaN). (stats.rs)
  * filter PASS iff: balance r=pop0/(pop0+pop1) in [0.2,0.8]; n>=10; d>=2;
    absmdiff>=0.2; AND permutation pval<=1e-3. (filter.rs)
  * permutation pval: hold matrix fixed, shuffle cluster labels nit=5000 times,
    pval = #(shuffled median-abs-diff >= observed)/nit; returns -1 (-> FAIL) if
    fewer than 8 reads. (pval.rs)

Two consequences the sweep figures must respect:
  * ASMS abstains (cannot PASS) below 10 reads by construction (n>=10, and the
    pval gate needs >=8 reads) -- a real design difference vs our DEPTH=2, not a bug.
  * The permutation gate multiplies per-region cost by ~nit stat recomputes; it is
    vectorized here (two matmuls), so it stays cheap.
"""

import numpy as np

PVAL_THR = 1e-3
PVAL_NIT = 5000
PVAL_MIN_READS = 8          # pval.rs: returns -1 if cl.len() < 8
FILTER_MIN_READS = 10       # filter.rs c2: n >= 10
FILTER_MIN_CPGS = 2         # filter.rs c2: d >= 2
FILTER_MIN_ABSMDIFF = 0.2   # filter.rs c2: absmdiff >= 0.2
FILTER_BALANCE = (0.2, 0.8) # filter.rs c1: r in [0.2, 0.8]
DEFAULT_MAXIT = 10          # cluster.rs: --maxit default


def _log_bernoulli_term(x, mu):
    """log(mu^x (1-mu)^(1-x)), matching prob.rs::log_bernoulli_term.

    No pseudocount. The (x==0,mu==0) and (x==1,mu==1) corners return log(1)=0
    exactly; everything else is x*ln(mu) + (1-x)*ln(1-mu) (may be -inf)."""
    x = np.asarray(x, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = x * np.log(mu) + (1.0 - x) * np.log1p(-mu)
    corner = ((x == 0.0) & (mu == 0.0)) | ((x == 1.0) & (mu == 1.0))
    return np.where(corner, 0.0, term)


def _logmvbern_rows(X, mu):
    """Per-read log multivariate-Bernoulli likelihood under cluster-mean vector mu.

    X: (n, d) with NaN for missing (Rust uses -1). mu: (d,). Missing CpGs are
    skipped (contribute 0), matching prob.rs::logmvbern."""
    obs = ~np.isnan(X)
    Xz = np.where(obs, X, 0.0)
    # ASMS uses unchecked f32: log(0)->-inf, NaN mu -> NaN, propagated silently.
    # Mirror that here (the BIC module set np.seterr(all='raise') at import).
    with np.errstate(all="ignore"):
        terms = _log_bernoulli_term(Xz, mu[None, :])  # (n, d)
        terms = np.where(obs, terms, 0.0)
        return terms.sum(axis=1)                       # (n,)


def _bernoulli_gamma(X, pi, mu):
    """E-step: responsibilities gamma (n, 2). cluster.rs::bernoulli_gamma.

    When a cluster has collapsed (pi[k]=0 or mu[k] NaN), log(0)/NaN propagate
    exactly as in ASMS's f32 code; the region then fails the balance filter."""
    n = X.shape[0]
    log_g = np.empty((n, 2))
    with np.errstate(all="ignore"):
        for k in range(2):
            log_g[:, k] = np.log(pi[k]) + _logmvbern_rows(X, mu[k])
        m = log_g.max(axis=1, keepdims=True)
        denom = m[:, 0] + np.log(np.exp(log_g - m).sum(axis=1))
        return np.exp(log_g - denom[:, None])


def _update_mu(X, gamma):
    """M-step: returns (mu (2,d), pop (2,)). cluster.rs::update_mu.

    mu[k,j] = sum_{i: obs} gamma[i,k]*X[i,j] / pop[k], pop[k]=sum_i gamma[i,k]
    (GLOBAL responsibility denominator, not per-position observed). An empty
    cluster (pop[k]=0) yields NaN mu, propagated as in ASMS's unchecked f32."""
    obs = ~np.isnan(X)
    Xz = np.where(obs, X, 0.0)
    with np.errstate(all="ignore"):
        pop = gamma.sum(axis=0)                        # (2,)
        num = gamma.T @ Xz                             # (2, d)
        mu = num / pop[:, None]
    return mu, pop


def cluster_matrix(X, seed, maxit=DEFAULT_MAXIT):
    """Run ASMS's K=2 Bernoulli EM. Returns (clusters (n,), mu (2,d), pop (2,)).

    Single random init (one restart), exactly `maxit` iterations, no convergence
    check. Final hard label: cluster 0 iff cprob0 > cprob1 (ties -> 1), matching
    cluster.rs's `if cprob[i][0] > cprob[i][1]`."""
    n, d = X.shape
    rng = np.random.RandomState(seed)
    labels = rng.randint(0, 2, size=n)                 # fastrand::u8(0..2) per read

    gamma0 = np.zeros((n, 2))
    gamma0[np.arange(n), labels] = 1.0
    pi = gamma0.sum(axis=0) / n                         # init_pi
    mu, _ = _update_mu(X, gamma0)                       # init_mu

    gamma = gamma0
    for _ in range(maxit):
        gamma = _bernoulli_gamma(X, pi, mu)
        mu, pop = _update_mu(X, gamma)
        pi = pop / n

    clusters = np.where(gamma[:, 0] > gamma[:, 1], 0, 1).astype(np.uint8)
    pop = np.array([(clusters == 0).sum(), (clusters == 1).sum()], dtype=float)
    return clusters, mu, pop


def _cluster_means(X, clusters):
    """Per-CpG methylated fraction for each cluster: (m0 (d,), m1 (d,)).

    m_c[j] = mean of observed 0/1 states among cluster-c reads at CpG j; NaN if
    no cluster-c read is observed at j. stats.rs m0/m1."""
    obs = ~np.isnan(X)
    Xz = np.where(obs, X, 0.0)
    out = []
    for c in (0, 1):
        sel = clusters == c
        cnt = obs[sel].sum(axis=0)
        meth = Xz[sel].sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            m = np.where(cnt > 0, meth / cnt, np.nan)
        out.append(m)
    return out[0], out[1]


def median_abs_diff(X, clusters):
    """absmdiff: median across CpGs of per-CpG |m1-m0|, dropping NaN CpGs.

    stats.rs::compute_stats median_abs_diff_notna_0_1. NaN if no valid CpG."""
    m0, m1 = _cluster_means(X, clusters)
    diff = np.abs(m1 - m0)
    diff = diff[~np.isnan(diff)]
    if diff.size == 0:
        return np.nan
    return float(np.median(diff))


def permutation_pval(X, clusters, nit=PVAL_NIT, rng=None):
    """Empirical p-value of the observed median-abs-diff under label shuffling.

    Returns (refmedian, pval). pval = #(shuffled median >= observed)/nit. Holds X
    fixed, shuffles cluster labels (preserving pop0/pop1). Returns pval=-1 if
    fewer than PVAL_MIN_READS reads. pval.rs::pval_of_matrix (nit=5000)."""
    n = X.shape[0]
    refmedian = median_abs_diff(X, clusters)
    if n < PVAL_MIN_READS:
        return refmedian, -1.0
    if rng is None:
        rng = np.random.RandomState(0)
    if np.isnan(refmedian):
        # Rust: `rndmedian >= NaN` is always false -> pval 0. (Region fails c2 anyway.)
        return refmedian, 0.0

    obs = ~np.isnan(X)
    Xz = np.where(obs, X, 0.0).astype(np.float64)
    obs_f = obs.astype(np.float64)

    # Monte-Carlo block: NaN-slice medians and 0/0 divisions are expected here;
    # suppress the global np.seterr(all='raise') set by EMfunctions_ASM.
    with np.errstate(all="ignore"):
        # S: (nit, n) one-hot for "read i in cluster 1" under each shuffle,
        # as permutations of the label vector (pop0/pop1 preserved).
        perm = np.argsort(rng.random((nit, n)), axis=1)
        lab1 = (clusters == 1).astype(np.float64)
        S = lab1[perm]                                 # (nit, n) in {0,1}
        notS = 1.0 - S

        cnt1 = S @ obs_f;   meth1 = S @ Xz             # (nit, d)
        cnt0 = notS @ obs_f; meth0 = notS @ Xz
        m1 = np.where(cnt1 > 0, meth1 / cnt1, np.nan)
        m0 = np.where(cnt0 > 0, meth0 / cnt0, np.nan)
        diff = np.abs(m1 - m0)                          # (nit, d)
        # median across CpGs ignoring NaN; all-NaN rows -> NaN (counts as < ref)
        rndmedians = np.nanmedian(diff, axis=1)
        ge = np.greater_equal(rndmedians, refmedian)   # NaN >= ref is False
    pval = float(np.count_nonzero(ge)) / nit
    return refmedian, pval


def score_asms(X, seed, nit=PVAL_NIT, rng=None, maxit=DEFAULT_MAXIT):
    """ASMS cluster+filter call on one read x CpG matrix X (NaN=missing).

    Returns a dict: call (0/1), absmdiff, pval, pop0, pop1, n, d, balance.
    call=1 iff the full ASMS filter PASSes (heuristics AND permutation gate)."""
    n, d = X.shape
    if n == 0 or d == 0:
        return dict(call=0, absmdiff=np.nan, pval=np.nan, pop0=0, pop1=0,
                    n=n, d=d, balance=np.nan)
    clusters, _mu, pop = cluster_matrix(X, seed, maxit=maxit)
    pop0, pop1 = int(pop[0]), int(pop[1])
    absmdiff, pval = permutation_pval(X, clusters, nit=nit, rng=rng)
    r = pop0 / (pop0 + pop1) if (pop0 + pop1) > 0 else np.nan
    c1 = (FILTER_BALANCE[0] <= r <= FILTER_BALANCE[1])
    c2 = (n >= FILTER_MIN_READS) and (d >= FILTER_MIN_CPGS) and \
         (not np.isnan(absmdiff)) and (absmdiff >= FILTER_MIN_ABSMDIFF)
    passed = bool(c1 and c2 and (0.0 <= pval <= PVAL_THR))
    return dict(call=int(passed), absmdiff=float(absmdiff), pval=float(pval),
                pop0=pop0, pop1=pop1, n=int(n), d=int(d),
                balance=float(r) if not np.isnan(r) else np.nan)
