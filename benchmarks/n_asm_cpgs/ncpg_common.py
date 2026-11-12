# ABOUTME: Shared generation core for the partial-ASM (n_asm CpGs) sensitivity sweep.
# ABOUTME: Builds regions where n of J=10 CpGs are ASM-differential and the rest SHARED; reused by all 3 method runners.
"""Partial-ASM region generation shared by the n_asm_cpgs sweep runners.

The question this analysis answers: how many of a region's J=10 CpGs must carry
the allele-specific (ASM) signal before each detector (EM/BIC, naive Fisher,
ASMS) calls the region ASM? We add a third truth axis, `n_asm` (0..10), to the
existing depth x noH grid. A region with `n_asm` ASM CpGs has:

  - CpGs 0..n_asm-1  : allele-differential (allele0 ~ hyper band, allele1 ~ hypo)
  - CpGs n_asm..J-1  : SHARED (both alleles equal, one uniform draw) — no signal

`n_asm = 0` reduces (distributionally) to the SHARED negative class; `n_asm = 10`
to the STATIC_ASM positive class. CpG *position* is irrelevant to every scorer
(BIC is per-CpG, Fisher pools all CpGs, ASMS takes the median across CpGs), so
"first n_asm" is equivalent to any size-n_asm subset.

This module touches NO production code. It reuses the production generator
(generate_synthetic), the production matrix build (EMfunctions_ASM.Reshape2Matrix),
and the sweep's tag-treatment + approach helpers (run_sweep) read-only, swapping
only the region-theta construction for the partial-ASM builder here.
"""

import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.join(_HERE, os.pardir)           # benchmarks/
_EM = os.path.join(_PARENT, os.pardir, "em")
_SYNTH = os.path.join(_EM, "synthetic")
for _p in (_PARENT, _EM, _SYNTH):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import EMfunctions_ASM as em          # noqa: E402
import generate_synthetic as gs       # noqa: E402
import run_sweep as rs                # noqa: E402  (read-only reuse of grid/treatment/approach helpers)

# Geometry + biology inherited from the production sweep (single source of truth).
J = rs.J                              # 10 CpGs per region == window size -> one verdict.
CHROM = rs.CHROM
SHARED_RANGE = rs.SHARED_RANGE        # (0.0, 1.0)
ASM_HYPER_RANGE = rs.ASM_HYPER_RANGE  # (0.65, 1.00)
ASM_HYPO_RANGE = rs.ASM_HYPO_RANGE    # (0.00, 0.35)

# Coarse grid for this analysis: depth 10..100 step 10, noH
# 0..1 step 0.10. n_asm sweeps 0..10 (0 = SHARED negative baseline for precision;
# 1..10 = the requested partial-ASM runs). Distinct seed base -> its own paired
# read space (BIC/Fisher/ASMS score byte-identical reads via stable_seed_ncpg).
COARSE_DEPTHS = list(range(10, 101, 10))
COARSE_NOHS = [round(i * 0.10, 2) for i in range(11)]
N_ASM_VALUES = list(range(0, 11))
DEFAULT_REPS = 100
SEED_BASE = 20260701


def coarse_grid():
    """Return (depths, nohs, default_reps) for the coarse n_asm sweep."""
    return list(COARSE_DEPTHS), list(COARSE_NOHS), DEFAULT_REPS


def n_asm_for_task(task_id):
    """Map a 1-based SGE array task id (1..11) to its n_asm value (0..10)."""
    if not (1 <= task_id <= len(N_ASM_VALUES)):
        raise SystemExit(
            f"--task-id {task_id} out of range 1..{len(N_ASM_VALUES)} (n_asm 0..10)"
        )
    return N_ASM_VALUES[task_id - 1]


def stable_seed_ncpg(seed_base, n_asm, depth, noh_idx, rep):
    """Deterministic 32-bit seed for one partial-ASM region.

    Same mixing form as run_sweep._stable_seed, with n_asm folded in where the
    class id used to live. Pure-Python integer arithmetic so the module-level
    np.seterr(all='raise') in EMfunctions_ASM cannot trip on overflow. All three
    method runners call THIS function, so for any (n_asm, depth, noH, rep) they
    score byte-identical reads (paired comparison).
    """
    s = (
        seed_base * 2654435761
        + n_asm * 40503
        + depth * 2246822519
        + noh_idx * 3266489917
        + rep * 668265263
    )
    return s & 0xFFFFFFFF


def partial_asm_theta(J, n_asm, rng,
                      shared_range=SHARED_RANGE,
                      asm_hyper_range=ASM_HYPER_RANGE,
                      asm_hypo_range=ASM_HYPO_RANGE):
    """True per-allele theta for a region with `n_asm` ASM CpGs (rest SHARED).

    Mirrors the distributions in generate_synthetic.generate_region_theta:
      - every CpG starts SHARED: both alleles = one uniform(shared_range) draw;
      - the first n_asm CpGs are overwritten with an allele difference at k=0
        (allele0 ~ uniform(asm_hyper_range), allele1 ~ uniform(asm_hypo_range));
      - k=1 is filled with free per-CpG values (unused by the sampler, pi[:,1]=0,
        and dead under the M1 diploid constraint) so the EM must learn "k=1
        unused" from pi, not from a degenerate theta.

    Returns theta of shape (J, 2, 2). n_asm=0 -> pure SHARED; n_asm=J -> full ASM.
    """
    theta = np.full((J, 2, 2), 0.5)
    # SHARED baseline at k=0 for every CpG (both alleles identical, no signal).
    vals = rng.uniform(shared_range[0], shared_range[1], J)
    theta[:, :, 0] = vals[:, None]
    # Overwrite the first n_asm CpGs with the allele-differential ASM signal.
    if n_asm > 0:
        theta[:n_asm, 0, 0] = rng.uniform(asm_hyper_range[0], asm_hyper_range[1], n_asm)
        theta[:n_asm, 1, 0] = rng.uniform(asm_hypo_range[0], asm_hypo_range[1], n_asm)
    # k=1: free per-CpG values (unused), matching the generator.
    theta[:, 0, 1] = rng.uniform(0.01, 0.99, J)
    theta[:, 1, 1] = rng.uniform(0.01, 0.99, J)
    return theta


def build_partial_region(depth, noh_frac, n_asm, rng):
    """Generate one partial-ASM region's reads -> (X, tags) after dropping all-NaN.

    Copy of run_sweep._build_region with the region-theta swapped for the
    partial-ASM builder; every other step (coords, read sampling, matrix build,
    all-NaN drop) is the production path, unchanged.
    """
    coords = gs.generate_coords(J, 1_000_000, rng)
    theta = partial_asm_theta(J, n_asm, rng)
    hap_probs = ((1 - noh_frac) / 2, (1 - noh_frac) / 2, noh_frac)
    reads = gs.generate_reads_for_region(
        coords, theta, depth, CHROM, rng, hap_probs=hap_probs,
    )
    df = pd.DataFrame(reads)
    coords_str = ",".join(df["startcoord"].astype(str))
    coords_array = np.array(sorted({int(float(c)) for c in coords_str.split(",")}))
    X = em.Reshape2Matrix(df, coords_array)
    tags = np.array(df["haplotype"])
    valid = ~np.isnan(X).all(axis=1)
    return X[valid], tags[valid]
