# ABOUTME: EM/BIC scorer over the partial-ASM (n_asm) x depth x noH coarse grid, 4 tag-treatments.
# ABOUTME: One SGE array task per n_asm value; reuses run_sweep._score (production BIC path) and ncpg_common generation.
"""EM/BIC arm of the n_asm-CpG sensitivity sweep.

For a fixed n_asm (selected by the SGE array task id), score every
(depth, noH, rep) partial-ASM region four ways (as_is / randomized / drop_noh /
all_noh) with the production BIC path (compute_null_BIC + _m1_best_of_n via
run_sweep._score). `n_asm = 0` regions are the SHARED negative baseline.

Output shard columns:
  n_asm, depth, noh_frac, rep, approach, total_reads,
  BICsinglecomp, BICmiddlecomp, BICresult, mean_theta_diff
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ncpg_common as nc              # noqa: E402
from ncpg_common import em, rs        # noqa: E402  (EMfunctions_ASM, run_sweep re-exported)

OUTPUT_COLUMNS = (
    "n_asm", "depth", "noh_frac", "rep", "approach", "total_reads",
    "BICsinglecomp", "BICmiddlecomp", "BICresult", "mean_theta_diff",
)


def run_cell(n_asm, depth, noh_frac, noh_idx, reps, seed_base, asm_restarts, tol):
    """Yield BIC rows for all reps x approaches of one (n_asm, depth, noH) cell."""
    rows = []
    for rep in range(reps):
        seed = nc.stable_seed_ncpg(seed_base, n_asm, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = nc.build_partial_region(depth, noh_frac, n_asm, rng)
        for approach in rs.APPROACHES:
            Xv, tv = rs._treatment(approach, X, tags, rng)
            n, bic0, bic1, result, mtd = rs._score(Xv, tv, asm_restarts, tol=tol)
            rows.append({
                "n_asm": n_asm,
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
    parser.add_argument("--task-id", type=int, required=True,
                        help="SGE array task id (1..11); selects n_asm = task_id-1 (0..10).")
    parser.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    parser.add_argument("--reps", type=int, default=None,
                        help="Regions per cell (default: ncpg_common.DEFAULT_REPS).")
    parser.add_argument("--asm-restarts", type=int, default=em.DEFAULT_RESTARTS)
    parser.add_argument("--em-tol", type=float, default=em.DEFAULT_EM_TOL,
                        help="ASM (M1) EM convergence tol (default %(default)g).")
    parser.add_argument("--seed-base", type=int, default=nc.SEED_BASE)
    args = parser.parse_args()

    depths, nohs, default_reps = nc.coarse_grid()
    reps = args.reps if args.reps is not None else default_reps
    n_asm = nc.n_asm_for_task(args.task_id)

    print(
        f"[bic] task={args.task_id}/{len(nc.N_ASM_VALUES)} n_asm={n_asm} "
        f"depth={len(depths)} noH={len(nohs)} reps={reps} "
        f"restarts={args.asm_restarts} em_tol={args.em_tol:g}"
    )

    all_rows = []
    for depth in depths:
        for noh_idx, noh_frac in enumerate(nohs):
            all_rows.extend(run_cell(
                n_asm, depth, noh_frac, noh_idx, reps,
                args.seed_base, args.asm_restarts, args.em_tol,
            ))

    out = pd.DataFrame(all_rows, columns=list(OUTPUT_COLUMNS))
    out = out.sort_values(["depth", "noh_frac", "rep", "approach"]).reset_index(drop=True)

    tmp = args.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.out)
    print(f"[bic] Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
