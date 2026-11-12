# ABOUTME: 2-way (Shared M0 vs ASM M1) EM-restart calibration sweep for the static ASM model.
# ABOUTME: Emits per-restart M1 BIC + the deterministic M0 BIC over a depth x noH grid so recall/FPR can be computed vs restart count.
"""EM-restart calibration for the static, 2-way (Shared vs ASM) model.

Reproduces the recall dip at high untagged-read fraction (noH ~ 0.7-0.8) that
the M1 EM shows when run from a single initialisation, and quantifies how many
random restarts remove it.

Cell generation and scoring are imported verbatim from `../run_sweep.py` (same
generator, ranges, seed function and seed_base). The regions are NOT the
benchmark grid's regions, though: noh_idx enumerates this sweep's own 0.05-step
noH list, so only noH = 0.0 shares seeds with the benchmark.

For each (truth_class in {SHARED, STATIC_ASM}, depth, noH) cell and each rep:
  * M0 (Shared) BIC once (closed form, deterministic).
  * M1 (ASM) BIC K times, restart r using the offset seed
    (_seed_from_X(X) + r*1_000_003) & 0xFFFFFFFF -- the exact scheme in
    EMfunctions_ASM._m1_best_of_n, so cumulative-min over restarts reproduces the
    production best-of-n. r=0 is bit-identical to the single-init default.

Output schema (long form, one row per (class, depth, noH, rep, restart_idx)):
  truth_class, depth, noh_frac, rep, approach, total_reads, restart_idx,
  BIC_M0, BIC_M1

The verdict at restart-count R is: ASM iff min(BIC_M1[0..R-1]) < BIC_M0.
Recall (STATIC_ASM) and FPR (SHARED) per (depth, noH, R) follow from that verdict.

SGE array: one task per depth (--task-id selects the depth from --depths).
"""

import argparse
import os
import sys
import time
from itertools import product
from multiprocessing import Pool

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_PP = os.path.normpath(os.path.join(_HERE, os.pardir))   # benchmarks/
if _PP not in sys.path:
    sys.path.insert(0, _PP)

import run_sweep as rs               # noqa: E402  (also puts em/ + synthetic/ on sys.path)

em = rs.em                            # EMfunctions_ASM
gs = rs.gs                            # generate_synthetic

OUTPUT_COLUMNS = (
    "truth_class", "depth", "noh_frac", "rep", "approach", "total_reads",
    "restart_idx", "BIC_M0", "BIC_M1",
)


def _per_restart_m1_bics(X, tags, K, tol):
    """M1 BIC for restart indices 0..K-1 using _m1_best_of_n's offset-seed scheme.

    Patches EMfunctions_ASM._seed_from_X (the name compute_M1_BIC resolves) so
    restart r draws a different theta init; restores it afterwards. The list's
    running minimum equals _m1_best_of_n(..., n_restarts=R) for every R<=K.
    """
    original = em._seed_from_X
    bics = []
    try:
        for r in range(K):
            em._seed_from_X = (
                lambda X, _r=r: (original(X) + _r * 1_000_003) & 0xFFFFFFFF
            )
            bic, _theta, _n_iters = em.compute_M1_BIC(X, tags, tol=tol)
            bics.append(float(bic))
    finally:
        em._seed_from_X = original
    return bics


def run_one_cell(cellargs):
    """One (truth_class, noH) cell at a fixed depth. Returns a long-form DataFrame."""
    (region_class, depth, noh_frac, noh_idx, reps, seed_base, K, tol, approach) = cellargs
    class_id = rs._CLASS_ID[region_class]
    rows = []
    for rep in range(reps):
        seed = rs._stable_seed(seed_base, class_id, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = rs._build_region(region_class, depth, noh_frac, rng)
        Xv, tv = rs._treatment(approach, X, tags, rng)

        keep = ~np.isnan(Xv).all(axis=1)
        Xs, ts = Xv[keep], tv[keep]
        n = int(Xs.shape[0])

        if n < em.DEPTH:
            # Degenerate window: mirror _score's sentinel (M0 wins). BIC_M0=0,
            # BIC_M1=1 at every restart -> min(1) < 0 is False -> Shared.
            bic0 = 0.0
            m1_bics = [1.0] * K
        else:
            bic0 = float(em.compute_null_BIC(Xs))
            m1_bics = _per_restart_m1_bics(Xs, ts, K, tol)

        for r, b in enumerate(m1_bics):
            rows.append({
                "truth_class": region_class,
                "depth": depth,
                "noh_frac": noh_frac,
                "rep": rep,
                "approach": approach,
                "total_reads": n,
                "restart_idx": r,
                "BIC_M0": bic0,
                "BIC_M1": b,
            })
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def _log_progress(done, total, t0):
    elapsed = time.time() - t0
    rate = done / elapsed if elapsed > 0 else 0
    eta = (total - done) / rate if rate > 0 else float("inf")
    print(f"  [{done}/{total}] elapsed={elapsed/60:.1f}min ETA={eta/60:.1f}min", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--depths", type=int, nargs="+", default=[30, 50, 100],
                    help="Depths to sweep. --task-id picks one for an SGE array.")
    ap.add_argument("--task-id", type=int, default=None,
                    help="1-based SGE array id; selects one depth from --depths.")
    ap.add_argument("--nohs", type=float, nargs="+", default=None,
                    help="noH fractions. Default: 0.00..1.00 step 0.05.")
    ap.add_argument("--reps", type=int, default=100, help="Regions per class per cell.")
    ap.add_argument("--K", type=int, default=10, help="EM restarts per region.")
    ap.add_argument("--classes", nargs="+", default=[gs.SHARED, gs.STATIC_ASM],
                    choices=[gs.SHARED, gs.STATIC_ASM])
    ap.add_argument("--approach", type=str, default="as_is", choices=list(rs.APPROACHES),
                    help="Tag treatment. Default as_is (Standard) -- where the noH dip lives.")
    ap.add_argument("--em-tol", type=float, default=em.DEFAULT_EM_TOL)
    ap.add_argument("--seed-base", type=int, default=20260603,
                    help="Match run_sweep default so regions align with the benchmark grid.")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    args = ap.parse_args()

    nohs = args.nohs if args.nohs is not None else [round(i * 0.05, 2) for i in range(21)]

    if args.task_id is not None:
        if not (1 <= args.task_id <= len(args.depths)):
            raise SystemExit(f"--task-id {args.task_id} out of range 1..{len(args.depths)}")
        depths = [args.depths[args.task_id - 1]]
    else:
        depths = args.depths

    cells = []
    for depth in depths:
        for noh_idx, noh_frac in enumerate(nohs):
            for region_class in args.classes:
                cells.append((
                    region_class, depth, noh_frac, noh_idx, args.reps,
                    args.seed_base, args.K, args.em_tol, args.approach,
                ))

    print(
        f"Restart sweep: depths={depths} nohs={len(nohs)} classes={args.classes} "
        f"reps={args.reps} K={args.K} approach={args.approach} "
        f"-> {len(cells)} cells, {args.workers} workers",
        flush=True,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    tmp_path = args.out + ".tmp"

    t0 = time.time()
    header_written = False
    completed = 0
    with open(tmp_path, "w") as f:
        if args.workers <= 1:
            for cell in cells:
                df = run_one_cell(cell)
                df.to_csv(f, sep="\t", index=False, header=not header_written)
                header_written = True
                f.flush()
                completed += 1
                _log_progress(completed, len(cells), t0)
        else:
            with Pool(processes=args.workers) as pool:
                for df in pool.imap_unordered(run_one_cell, cells):
                    df.to_csv(f, sep="\t", index=False, header=not header_written)
                    header_written = True
                    f.flush()
                    completed += 1
                    _log_progress(completed, len(cells), t0)

    os.replace(tmp_path, args.out)
    print(f"Done. {completed}/{len(cells)} cells in {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
