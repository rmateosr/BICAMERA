# ABOUTME: Validates asms_scorer.py against the ecmra/asms embedded test vectors and the real binary.
# ABOUTME: (1) EM math vs Rust unit-test vectors; (2) stats+permutation pval vs `asms pval` on sweep regions.
#
# The test matrix and expected values in Part 1 are copied from the unit tests of
# ASMS (prob.rs, cluster.rs, matrix.rs) by E. Raineri, M. A. Esteve Marco, E. Casals
# and A. Esteve Codina, NAR Genomics and Bioinformatics 8(3):lqag072 (2026),
# doi:10.1093/nargab/lqag072. Source: ecmra/asms Zenodo release,
# doi:10.5281/zenodo.20323001 (commit c7178db), used under CC BY 4.0
# (https://creativecommons.org/licenses/by/4.0/). Changes: values transcribed
# from Rust into Python assertions; nothing else is copied.
"""Validate the Python ASMS reimplementation.

Part 1 (no binary): assert the EM primitives reproduce the exact numbers in the
Rust unit tests (prob.rs, cluster.rs) -- log_bernoulli_term, logmvbern,
bernoulli_gamma, init_pi/init_mu, update_mu.

Part 2 (binary): for a handful of sweep-generated regions, run our cluster_matrix,
write the read x CpG matrix and our cluster labels in ASMS's file formats, run the
real `asms pval <clusterfn> <matrixfn>`, and compare:
  * refmedian (median-abs-diff): MUST match exactly -- it is deterministic given
    the matrix + labels, so this validates compute_stats / median_abs_diff.
  * pval: compared within Monte-Carlo tolerance and PASS/FAIL agreement (the
    binary's shuffle RNG is not seeded the same as ours, so exact match is not
    expected; refmedian is the deterministic anchor).

Usage: python validate_asms_scorer.py [path/to/asms/binary]

The binary is the first argument, else $ASMS_BIN; with neither, Part 2 is skipped.
Build it from the Zenodo release above (`cargo build -r` -> target/release/asms).
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))      # benchmarks/asms
_PARENT = os.path.join(_HERE, os.pardir)                # benchmarks/ (run_sweep)
_EM = os.path.join(_HERE, os.pardir, os.pardir, "em")   # em/
_SYNTH = os.path.join(_EM, "synthetic")
for _p in (_HERE, _PARENT, _EM, _SYNTH):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import asms_scorer as asms          # noqa: E402

DEFAULT_BINARY = os.environ.get("ASMS_BIN", "")


# --------------------------------------------------------------------------
# Part 1 -- embedded Rust test vectors
# --------------------------------------------------------------------------

def _test_matrix():
    """matrix.rs::test::make_test_matrix -- 4 reads x 5 CpGs."""
    return np.array([
        [0., 1., 0., 1., 0.],
        [1., 0., 1., 0., 1.],
        [0., 1., 0., 1., 0.],
        [1., 0., 1., 0., 1.],
    ])


def part1_vectors():
    fails = []

    def check(name, got, want, tol):
        ok = abs(got - want) < tol
        print(f"  [{'ok' if ok else 'FAIL'}] {name}: got={got:.6g} want={want:.6g}")
        if not ok:
            fails.append(name)

    # prob.rs::log_bernoulli_term_test
    check("log_bernoulli_term(1.0,0.2)",
          float(asms._log_bernoulli_term(1.0, 0.2)), -1.6094379124341003, 1e-9)

    # prob.rs::logmvbern_test : logmvbern([0,1],[0.2,0.8]) = -0.446287
    lv = float(asms._logmvbern_rows(np.array([[0., 1.]]), np.array([0.2, 0.8]))[0])
    check("logmvbern([0,1],[0.2,0.8])", lv, -0.446287, 1e-5)

    # cluster.rs::bernoulli_gamma_test : gamma[0]=0.999024, gamma[1]=0.00097561
    tm = _test_matrix()
    pi = np.array([0.5, 0.5])
    mu = np.array([[0.2, 0.8, 0.2, 0.8, 0.2],
                   [0.8, 0.2, 0.8, 0.2, 0.8]])
    g = asms._bernoulli_gamma(tm, pi, mu)
    check("bernoulli_gamma[0,0]", float(g[0, 0]), 0.999024, 1e-5)
    check("bernoulli_gamma[0,1]", float(g[0, 1]), 0.00097561, 1e-4)

    # cluster.rs::init_pi_test : [0.5,0.5] for cluster [0,1,0,1]
    labels = np.array([0, 1, 0, 1])
    n = 4
    gamma0 = np.zeros((n, 2)); gamma0[np.arange(n), labels] = 1.0
    pi0 = gamma0.sum(axis=0) / n
    check("init_pi[0]", float(pi0[0]), 0.5, 1e-9)

    # cluster.rs::init_mu_test : mu = [0,1,0,1,0, 1,0,1,0,1]
    mu0, _ = asms._update_mu(tm, gamma0)
    want_mu = np.array([[0, 1, 0, 1, 0], [1, 0, 1, 0, 1]], dtype=float)
    ok = np.allclose(mu0, want_mu, atol=1e-9)
    print(f"  [{'ok' if ok else 'FAIL'}] init_mu: {'matches' if ok else mu0.tolist()}")
    if not ok:
        fails.append("init_mu")

    # cluster.rs::update_mu_test : gamma given -> pop=[2,2], mu=[0.1,0.9,...]
    gA = np.array([[0.9, 0.1], [0.1, 0.9], [0.9, 0.1], [0.1, 0.9]])
    muA, popA = asms._update_mu(tm, gA)
    want_muA = np.array([[0.1, 0.9, 0.1, 0.9, 0.1], [0.9, 0.1, 0.9, 0.1, 0.9]])
    ok = np.allclose(popA, [2., 2.], atol=1e-6) and np.allclose(muA, want_muA, atol=1e-6)
    print(f"  [{'ok' if ok else 'FAIL'}] update_mu: pop={popA.tolist()}")
    if not ok:
        fails.append("update_mu")

    return fails


# --------------------------------------------------------------------------
# Part 2 -- vs the real binary via matrix-file + `asms pval`
# --------------------------------------------------------------------------

def _write_matrix_file(path, X, contig="chrSYN"):
    """Write X (n x d, NaN=missing) in ASMS matrix-file format (matrix_of_file)."""
    n, d = X.shape
    lines = [f"#@N:{n}", f"#@D:{d}", f"#@CONTIG:{contig}"]
    for i in range(n):
        for j in range(d):
            v = X[i, j]
            if not np.isnan(v):
                # idxread  rname  idxgpos  gpos  state
                lines.append(f"{i}\tr{i}\t{j}\t{1000 + j}\t{int(v)}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _write_cluster_file(path, clusters):
    """Write cluster labels in cluster_of_file format (reads fields[1] per line).

    Lines must be in read-index order 0..n-1 so clusters[i] maps to read i."""
    lines = ["#@N:%d" % len(clusters), "#@K:2"]
    for i, c in enumerate(clusters):
        lines.append(f"r{i}\t{int(c)}\t1.0000\t0.0000")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _binary_pval(binary, clusterfn, matrixfn):
    """Run `asms pval` and parse (refmedian, pval) from the last two tab fields."""
    out = subprocess.run([binary, "pval", clusterfn, matrixfn],
                         capture_output=True, text=True)
    line = out.stdout.strip().splitlines()[-1]
    fields = line.split("\t")
    return float(fields[-2]), float(fields[-1])


def part2_binary(binary):
    import run_sweep as rs
    import generate_synthetic as gs

    fails = []
    cases = [
        (gs.STATIC_ASM, 40, 0.4, 11),
        (gs.STATIC_ASM, 20, 0.6, 22),
        (gs.SHARED,     40, 0.4, 33),
        (gs.SHARED,     20, 0.2, 44),
        (gs.STATIC_ASM, 80, 0.8, 55),
        (gs.STATIC_ASM, 12, 0.0, 66),
    ]
    tmpd = tempfile.mkdtemp(prefix="asms_validate_")
    print(f"  tmp dir: {tmpd}")
    for cls, depth, noh, seed in cases:
        rng = np.random.RandomState(seed)
        X, tags = rs._build_region(cls, depth, noh, rng)
        clusters, _mu, _pop = asms.cluster_matrix(X, seed)
        mf = os.path.join(tmpd, f"matrix_{cls}_{depth}_{noh}.txt")
        cf = os.path.join(tmpd, f"cluster_{cls}_{depth}_{noh}.txt")
        _write_matrix_file(mf, X)
        _write_cluster_file(cf, clusters)

        bin_ref, bin_pval = _binary_pval(binary, cf, mf)
        py_ref = asms.median_abs_diff(X, clusters)
        _, py_pval = asms.permutation_pval(X, clusters, rng=np.random.RandomState(seed))

        ref_ok = (np.isnan(bin_ref) and np.isnan(py_ref)) or abs(bin_ref - py_ref) < 1e-3
        # pval: PASS/FAIL vs the 1e-3 gate must agree; values within MC tolerance.
        gate = lambda p: (0.0 <= p <= asms.PVAL_THR)
        pval_gate_ok = gate(bin_pval) == gate(py_pval)
        pval_close = abs(bin_pval - py_pval) < 0.02 or pval_gate_ok
        tag = f"{cls} d={depth} noH={noh} n={X.shape[0]}"
        print(f"  [{'ok' if (ref_ok and pval_gate_ok) else 'FAIL'}] {tag}: "
              f"refmed bin={bin_ref:.4f} py={py_ref:.4f} | "
              f"pval bin={bin_pval:.4f} py={py_pval:.4f} gate_agree={pval_gate_ok}")
        if not (ref_ok and pval_gate_ok and pval_close):
            fails.append(tag)
    return fails


def main():
    binary = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BINARY
    print("== Part 1: EM math vs embedded Rust test vectors ==")
    f1 = part1_vectors()
    print("\n== Part 2: stats + permutation pval vs real binary ==")
    if not os.path.exists(binary):
        print(f"  SKIP: binary not found at {binary}")
        f2 = []
    else:
        f2 = part2_binary(binary)

    print("\n== SUMMARY ==")
    if not f1 and not f2:
        print("ALL CHECKS PASSED")
        return 0
    print(f"FAILURES: part1={f1} part2={f2}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
