# ABOUTME: noH x depth sweep harness for the static ASM model, scoring 4 tag-treatments.
# ABOUTME: One generated dataset per (class, depth, noH, rep); emits per-region BIC rows for P/R.
"""Sweep harness: precision/recall of the static ASM BIC under 4 tag-treatments.

For each grid cell (depth, noH) and each truth class (SHARED negative,
STATIC_ASM positive), we generate `reps` independent regions. Each region's
reads are scored four ways on the SAME underlying reads:

  as_is       -- H1/H2/noH tags as generated.
  randomized  -- the tag vector permuted among reads (counts kept, link broken).
  drop_noh    -- noH reads discarded; only H1/H2 scored.
  all_noh     -- H1/H2 tags stripped; every read scored as noH (fully unphased).

Scoring reuses the production code path exactly (EMfunctions_ASM):
compute_null_BIC (Shared/M0) and _m1_best_of_n (ASM/M1, 5 restarts), with
J = windowsize = 10 so each region is a single window -> one BIC verdict.
BICresult = 1 iff ASM (M1) beats Shared (M0).

SGE sharding: one array task per depth value (--task-id = SGE_TASK_ID).
Output: one TSV shard with columns
  truth_class, depth, noh_frac, rep, approach, total_reads,
  BICsinglecomp, BICmiddlecomp, BICresult, mean_theta_diff
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_EM = os.path.join(_HERE, os.pardir, "em")
_SYNTH = os.path.join(_EM, "synthetic")
for _p in (_EM, _SYNTH):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import EMfunctions_ASM as em          # noqa: E402
import generate_synthetic as gs       # noqa: E402

# Paper methylation ranges (Nature s41467-025-57433-1): ASM alleles at
# 0-35% and 65-100%. Negatives (SHARED) span the full range, equal on both
# alleles. Exposed as module constants so future runs can sweep other ranges.
ASM_HYPER_RANGE = (0.65, 1.00)
ASM_HYPO_RANGE = (0.00, 0.35)
SHARED_RANGE = (0.00, 1.00)

J = 10                 # CpGs per region == window size -> one verdict per region.
WINDOWSIZE = 10
CHROM = "chrSYN"

APPROACHES = ("as_is", "randomized", "drop_noh", "all_noh")

_CLASS_ID = {gs.SHARED: 0, gs.STATIC_ASM: 1}

OUTPUT_COLUMNS = (
    "truth_class", "depth", "noh_frac", "rep", "approach", "total_reads",
    "BICsinglecomp", "BICmiddlecomp", "BICresult", "mean_theta_diff",
)


GRIDS = ("smoke", "fine", "paper")


def grid(name):
    """Return (depths, noh_fracs, default_reps) for a named grid.

    paper is the manuscript's grid: depth 10-100 step 10 x noH 0.0-1.0 step 0.1,
    an exact subset of fine. Use noh_cells() to iterate it, not enumerate().
    """
    if name == "smoke":
        depths = [5, 10, 20, 40, 80]
        nohs = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        return depths, nohs, 50
    if name == "fine":
        depths = list(range(1, 101))
        nohs = [round(i * 0.01, 2) for i in range(101)]
        return depths, nohs, 100
    if name == "paper":
        depths = list(range(10, 101, 10))
        nohs = [round(i * 0.10, 2) for i in range(11)]
        return depths, nohs, 100
    raise ValueError(f"unknown grid: {name!r} (expected one of {GRIDS})")


def noh_cells(name):
    """Return the (noh_idx, noh_frac) pairs a grid seeds its regions with.

    noh_idx enters _stable_seed. On paper it is the noH's index on the fine
    0.01-step grid, so each paper cell draws the same regions as the matching
    fine cell; re-enumerating 0..10 would draw different ones. smoke and fine
    use their own enumeration (for fine the two are the same).
    """
    _, nohs, _ = grid(name)
    if name == "paper":
        return [(round(f * 100), f) for f in nohs]
    return list(enumerate(nohs))


def _stable_seed(seed_base, class_id, depth, noh_idx, rep):
    """Deterministic 32-bit seed for one region.

    Pure-Python integer arithmetic (no numpy) so module-level
    np.seterr(all='raise') in EMfunctions_ASM cannot trip on overflow.
    """
    s = (
        seed_base * 2654435761
        + class_id * 40503
        + depth * 2246822519
        + noh_idx * 3266489917
        + rep * 668265263
    )
    return s & 0xFFFFFFFF


def _build_region(region_class, depth, noh_frac, rng):
    """Generate one region's reads -> (X, tags) after dropping all-NaN reads."""
    coords = gs.generate_coords(J, 1_000_000, rng)
    theta = gs.generate_region_theta(
        region_class, J, rng,
        shared_range=SHARED_RANGE,
        asm_hyper_range=ASM_HYPER_RANGE,
        asm_hypo_range=ASM_HYPO_RANGE,
    )
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


def _treatment(approach, X, tags, rng):
    """Derive the (X, tags) variant for a tag-treatment from the base reads."""
    if approach == "as_is":
        return X, tags
    if approach == "randomized":
        return X, rng.permutation(tags)
    if approach == "drop_noh":
        keep = tags != "noH"
        return X[keep], tags[keep]
    if approach == "all_noh":
        return X, np.array(["noH"] * len(tags))
    raise ValueError(f"unknown approach: {approach!r}")


def _score(X, tags, asm_restarts, tol=em.DEFAULT_EM_TOL):
    """Score one (X, tags) -> (n_reads, BICsinglecomp, BICmiddlecomp, result, mtd).

    Mirrors EMBIC_bin_path_ASM's per-window rule: reads below DEPTH get the
    production sentinel (BICsinglecomp=0, BICmiddlecomp=1 -> Shared wins).
    """
    keep = ~np.isnan(X).all(axis=1)
    Xs = X[keep]
    ts = tags[keep]
    n = int(Xs.shape[0])
    if n < em.DEPTH:
        return n, 0.0, 1.0, 0, float("nan")
    bic_null = em.compute_null_BIC(Xs)
    bic_m1, theta, _ = em._m1_best_of_n(Xs, ts, asm_restarts, tol=tol)
    result = int(bic_m1 < bic_null)
    mtd = float(np.mean(np.abs(theta[:, 1, 0] - theta[:, 0, 0])))
    return n, float(bic_null), float(bic_m1), result, mtd


def run_cell(region_class, depth, noh_frac, noh_idx, reps, seed_base, asm_restarts,
             tol=em.DEFAULT_EM_TOL):
    """Yield rows for all reps x approaches of one (class, depth, noH) cell."""
    rows = []
    class_id = _CLASS_ID[region_class]
    for rep in range(reps):
        seed = _stable_seed(seed_base, class_id, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = _build_region(region_class, depth, noh_frac, rng)
        for approach in APPROACHES:
            Xv, tv = _treatment(approach, X, tags, rng)
            n, bic0, bic1, result, mtd = _score(Xv, tv, asm_restarts, tol=tol)
            rows.append({
                "truth_class": region_class,
                "depth": depth,
                "noh_frac": noh_frac,
                "rep": rep,
                "approach": approach,
                "total_reads": n,
                "BICsinglecomp": bic0,
                "BICmiddlecomp": bic1,
                "BICresult": result,
                "mean_theta_diff": mtd,
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", choices=GRIDS, required=True)
    parser.add_argument("--task-id", type=int, required=True,
                        help="SGE array task id (1-based); selects one depth.")
    parser.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    parser.add_argument("--reps", type=int, default=None,
                        help="Regions per class per cell (default: grid-specific).")
    parser.add_argument("--asm-restarts", type=int, default=em.DEFAULT_RESTARTS)
    parser.add_argument("--em-tol", type=float, default=em.DEFAULT_EM_TOL,
                        help="ASM (M1) EM convergence tol (default %(default)g).")
    parser.add_argument("--seed-base", type=int, default=20260603)
    args = parser.parse_args()

    depths, nohs, default_reps = grid(args.grid)
    reps = args.reps if args.reps is not None else default_reps

    if not (1 <= args.task_id <= len(depths)):
        raise SystemExit(
            f"--task-id {args.task_id} out of range 1..{len(depths)} for grid {args.grid!r}"
        )
    depth = depths[args.task_id - 1]

    print(
        f"grid={args.grid} task={args.task_id}/{len(depths)} depth={depth} "
        f"noH={len(nohs)} reps/class={reps} restarts={args.asm_restarts} "
        f"em_tol={args.em_tol:g}"
    )

    all_rows = []
    for noh_idx, noh_frac in noh_cells(args.grid):
        for region_class in (gs.SHARED, gs.STATIC_ASM):
            all_rows.extend(run_cell(
                region_class, depth, noh_frac, noh_idx, reps,
                args.seed_base, args.asm_restarts, tol=args.em_tol,
            ))

    out = pd.DataFrame(all_rows, columns=list(OUTPUT_COLUMNS))
    out = out.sort_values(["truth_class", "noh_frac", "rep", "approach"]).reset_index(drop=True)

    tmp = args.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.out)
    print(f"Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
