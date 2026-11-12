# ABOUTME: ASMS (Raineri et al. 2026) scorer over the partial-ASM (n_asm) x depth x noH coarse grid, 4 tag-treatments.
# ABOUTME: Reuses asms/asms_scorer.score_asms on the SAME reads as the BIC/Fisher arms (paired via stable_seed_ncpg).
"""ASMS arm of the n_asm-CpG sensitivity sweep.

Scores each partial-ASM region with the faithful ASMS cluster+filter
reimplementation (K=2 Bernoulli EM + balance / n>=10 / >=2 CpGs / median-abs-diff
>=0.2 / 5000-shuffle permutation p<=1e-3). ASMS is tag-agnostic, so as_is,
randomized, and all_noh see identical reads (identical call); only drop_noh
changes the read set. Its median-abs-diff gate is expected to require a MAJORITY
of CpGs to be ASM before it fires — the crux of this sweep.

Output shard columns:
  n_asm, depth, noh_frac, rep, approach, asms_n_reads, asms_call, asms_absmdiff, asms_pval
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

_ASMS = os.path.join(nc._PARENT, "asms")
if _ASMS not in sys.path:
    sys.path.insert(0, _ASMS)
import asms_scorer as asms            # noqa: E402

OUTPUT_COLUMNS = (
    "n_asm", "depth", "noh_frac", "rep", "approach",
    "asms_n_reads", "asms_call", "asms_absmdiff", "asms_pval",
)


def _asms_score(X, seed, nit):
    """ASMS cluster+filter on the treatment's reads (tag-agnostic).

    Mirrors run_asms_sweep._asms_score: EM init from `seed`, permutation shuffles
    from seed+1; drop reads with no observed CpG first (as the BIC path does)."""
    keep = ~np.isnan(X).all(axis=1)
    Xs = X[keep]
    if Xs.shape[0] == 0:
        return 0, 0, float("nan"), float("nan")
    res = asms.score_asms(Xs, seed=seed, nit=nit,
                          rng=np.random.RandomState((seed + 1) & 0xFFFFFFFF))
    return res["n"], res["call"], res["absmdiff"], res["pval"]


def run_cell(n_asm, depth, noh_frac, noh_idx, reps, seed_base, nit):
    """Yield ASMS rows for all reps x approaches of one (n_asm, depth, noH) cell."""
    rows = []
    for rep in range(reps):
        seed = nc.stable_seed_ncpg(seed_base, n_asm, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = nc.build_partial_region(depth, noh_frac, n_asm, rng)
        for approach in rs.APPROACHES:
            Xv, tv = rs._treatment(approach, X, tags, rng)
            n, call, absmdiff, pval = _asms_score(Xv, seed, nit)
            rows.append({
                "n_asm": n_asm,
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
    parser.add_argument("--task-id", type=int, required=True,
                        help="SGE array task id (1..11); selects n_asm = task_id-1 (0..10).")
    parser.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    parser.add_argument("--reps", type=int, default=None,
                        help="Regions per cell (default: ncpg_common.DEFAULT_REPS).")
    parser.add_argument("--asms-nit", type=int, default=asms.PVAL_NIT,
                        help="ASMS permutation shuffles per region (default 5000).")
    parser.add_argument("--seed-base", type=int, default=nc.SEED_BASE)
    args = parser.parse_args()

    depths, nohs, default_reps = nc.coarse_grid()
    reps = args.reps if args.reps is not None else default_reps
    n_asm = nc.n_asm_for_task(args.task_id)

    print(
        f"[asms] task={args.task_id}/{len(nc.N_ASM_VALUES)} n_asm={n_asm} "
        f"depth={len(depths)} noH={len(nohs)} reps={reps} nit={args.asms_nit} "
        f"seed_base={args.seed_base}"
    )

    all_rows = []
    for depth in depths:
        for noh_idx, noh_frac in enumerate(nohs):
            all_rows.extend(run_cell(
                n_asm, depth, noh_frac, noh_idx, reps, args.seed_base, args.asms_nit,
            ))

    out = pd.DataFrame(all_rows, columns=list(OUTPUT_COLUMNS))
    out = out.sort_values(["depth", "noh_frac", "rep", "approach"]).reset_index(drop=True)

    tmp = args.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.out)
    print(f"[asms] Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
