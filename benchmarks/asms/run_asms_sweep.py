# ABOUTME: ASMS (Raineri et al. 2026) ASM detector swept over depth x noH x 4 tag-treatments.
# ABOUTME: Fully separate from the EM/BIC harness; reuses its generation+seeds for a paired comparison.
"""ASMS baseline for the noH x depth ASM sweep — same framework as the Fisher run.

Per region, score the reads with a faithful Python reimplementation of ASMS's
`cluster`+`filter` (see asms_scorer.py): K=2 multivariate-Bernoulli EM (10 iters,
free pi, single init) then the ASMS filter (balance 20-80%, n>=10, >=2 CpGs,
median-abs-diff>=0.2, AND a 5000-shuffle permutation p<=1e-3). The call is binary
(`asms_call`); we also store the diagnostics (absmdiff, pval, n).

ASMS is **tag-agnostic**: its clustering ignores the haplotype tags. Under the four
tag-treatments the only thing that changes is the read SET, so:
  as_is == randomized == all_noh   (identical reads -> identical ASMS call)
  drop_noh                          (noH reads removed -> fewer reads -> may differ)
This is the point of the comparison: it shows ASMS gains nothing from tags, unlike
Fisher (which collapses to no-call under all_noh) and our BIC (which exploits them).

This pipeline is independent of run_sweep.py: it imports run_sweep ONLY to reuse
its generation, treatment, seeding, and grid helpers, so it scores the SAME reads
(same seed_base) as the BIC run — a paired comparison — without modifying it.
Output: one TSV shard with columns
  truth_class, depth, noh_frac, rep, approach, asms_n_reads, asms_call,
  asms_absmdiff, asms_pval
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.join(_HERE, os.pardir)          # benchmarks/
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import run_sweep as rs        # noqa: E402  (read-only reuse of generation/treatment/seeds)
import asms_scorer as asms    # noqa: E402

# MUST match run_sweep.py's --seed-base default so the generated reads are
# byte-identical to the BIC run (paired comparison). If you change one, change both.
SEED_BASE = 20260603

OUTPUT_COLUMNS = (
    "truth_class", "depth", "noh_frac", "rep", "approach",
    "asms_n_reads", "asms_call", "asms_absmdiff", "asms_pval",
)


def _asms_score(X, seed, nit):
    """ASMS cluster+filter on the treatment's reads. Returns (n, call, absmdiff, pval).

    Tag-agnostic (ignores `tags`); seed drives the EM init, seed+1 the permutation
    shuffles. Reads with no observed CpG are dropped first (as the BIC path does)."""
    keep = ~np.isnan(X).all(axis=1)
    Xs = X[keep]
    if Xs.shape[0] == 0:
        return 0, 0, float("nan"), float("nan")
    res = asms.score_asms(Xs, seed=seed, nit=nit,
                          rng=np.random.RandomState((seed + 1) & 0xFFFFFFFF))
    return res["n"], res["call"], res["absmdiff"], res["pval"]


def run_cell(region_class, depth, noh_frac, noh_idx, reps, seed_base, nit):
    """Yield ASMS rows for all reps x approaches of one (class, depth, noH) cell.

    Mirrors run_sweep.run_cell's generation/treatment sequence exactly (same
    seed -> same reads), swapping the EM/BIC scorer for the ASMS scorer."""
    rows = []
    class_id = rs._CLASS_ID[region_class]
    for rep in range(reps):
        seed = rs._stable_seed(seed_base, class_id, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = rs._build_region(region_class, depth, noh_frac, rng)
        for approach in rs.APPROACHES:
            Xv, tv = rs._treatment(approach, X, tags, rng)
            n, call, absmdiff, pval = _asms_score(Xv, seed, nit)
            rows.append({
                "truth_class": region_class,
                "depth": depth,
                "noh_frac": noh_frac,
                "rep": rep,
                "approach": approach,
                "asms_n_reads": n,
                "asms_call": call,
                "asms_absmdiff": absmdiff,
                "asms_pval": pval,
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", choices=rs.GRIDS, required=True)
    parser.add_argument("--task-id", type=int, required=True,
                        help="SGE array task id (1-based); selects one depth.")
    parser.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    parser.add_argument("--reps", type=int, default=None,
                        help="Regions per class per cell (default: grid-specific).")
    parser.add_argument("--asms-nit", type=int, default=asms.PVAL_NIT,
                        help="ASMS permutation shuffles per region (default 5000).")
    parser.add_argument("--seed-base", type=int, default=SEED_BASE,
                        help="Must equal run_sweep.py's --seed-base for paired reads.")
    args = parser.parse_args()

    depths, nohs, default_reps = rs.grid(args.grid)
    reps = args.reps if args.reps is not None else default_reps

    if not (1 <= args.task_id <= len(depths)):
        raise SystemExit(
            f"--task-id {args.task_id} out of range 1..{len(depths)} for grid {args.grid!r}"
        )
    depth = depths[args.task_id - 1]

    print(
        f"[asms] grid={args.grid} task={args.task_id}/{len(depths)} depth={depth} "
        f"noH={len(nohs)} reps/class={reps} nit={args.asms_nit} seed_base={args.seed_base}"
    )

    all_rows = []
    for noh_idx, noh_frac in rs.noh_cells(args.grid):
        for region_class in (rs.gs.SHARED, rs.gs.STATIC_ASM):
            all_rows.extend(run_cell(
                region_class, depth, noh_frac, noh_idx, reps, args.seed_base, args.asms_nit,
            ))

    out = pd.DataFrame(all_rows, columns=list(OUTPUT_COLUMNS))
    out = out.sort_values(["truth_class", "noh_frac", "rep", "approach"]).reset_index(drop=True)

    tmp = args.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.out)
    print(f"[asms] Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
