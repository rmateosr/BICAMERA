# ABOUTME: Naive Fisher-exact scorer over the partial-ASM (n_asm) x depth x noH coarse grid, 4 tag-treatments.
# ABOUTME: Reuses fisher/run_fisher_sweep._fisher_score on the SAME reads as the BIC/ASMS arms (paired via stable_seed_ncpg).
"""Fisher-exact arm of the n_asm-CpG sensitivity sweep.

Pools each region's H1/H2 CpG calls into one 2x2 (allele x methylation) table and
stores the raw two-sided Fisher p-value (thresholded at plot time). noH reads are
excluded — the naive baseline's defining limitation. Scores the identical reads
the BIC arm sees (same ncpg_common generation + stable_seed_ncpg).

Output shard columns:
  n_asm, depth, noh_frac, rep, approach, fisher_n_reads, fisher_p
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
from ncpg_common import rs            # noqa: E402  (run_sweep, for APPROACHES/_treatment)

_FISHER = os.path.join(nc._PARENT, "fisher")
if _FISHER not in sys.path:
    sys.path.insert(0, _FISHER)
from run_fisher_sweep import _fisher_score  # noqa: E402

OUTPUT_COLUMNS = (
    "n_asm", "depth", "noh_frac", "rep", "approach", "fisher_n_reads", "fisher_p",
)


def run_cell(n_asm, depth, noh_frac, noh_idx, reps, seed_base):
    """Yield Fisher rows for all reps x approaches of one (n_asm, depth, noH) cell."""
    rows = []
    for rep in range(reps):
        seed = nc.stable_seed_ncpg(seed_base, n_asm, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = nc.build_partial_region(depth, noh_frac, n_asm, rng)
        for approach in rs.APPROACHES:
            Xv, tv = rs._treatment(approach, X, tags, rng)
            n_used, p = _fisher_score(Xv, tv)
            rows.append({
                "n_asm": n_asm,
                "depth": depth,
                "noh_frac": noh_frac,
                "rep": rep,
                "approach": approach,
                "fisher_n_reads": n_used,
                "fisher_p": p,
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-id", type=int, required=True,
                        help="SGE array task id (1..11); selects n_asm = task_id-1 (0..10).")
    parser.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    parser.add_argument("--reps", type=int, default=None,
                        help="Regions per cell (default: ncpg_common.DEFAULT_REPS).")
    parser.add_argument("--seed-base", type=int, default=nc.SEED_BASE)
    args = parser.parse_args()

    depths, nohs, default_reps = nc.coarse_grid()
    reps = args.reps if args.reps is not None else default_reps
    n_asm = nc.n_asm_for_task(args.task_id)

    print(
        f"[fisher] task={args.task_id}/{len(nc.N_ASM_VALUES)} n_asm={n_asm} "
        f"depth={len(depths)} noH={len(nohs)} reps={reps} seed_base={args.seed_base}"
    )

    all_rows = []
    for depth in depths:
        for noh_idx, noh_frac in enumerate(nohs):
            all_rows.extend(run_cell(n_asm, depth, noh_frac, noh_idx, reps, args.seed_base))

    out = pd.DataFrame(all_rows, columns=list(OUTPUT_COLUMNS))
    out = out.sort_values(["depth", "noh_frac", "rep", "approach"]).reset_index(drop=True)

    tmp = args.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.out)
    print(f"[fisher] Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
