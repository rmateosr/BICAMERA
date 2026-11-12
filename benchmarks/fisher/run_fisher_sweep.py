# ABOUTME: Naive Fisher-exact ASM detector swept over depth x noH x 4 tag-treatments.
# ABOUTME: Fully separate from the EM/BIC harness; reuses its generation+seeds for a paired comparison.
"""Fisher-exact baseline for the noH x depth ASM sweep.

A deliberately naive comparison to the EM/BIC method: per region, pool all CpG
methylation calls into a 2x2 (allele x methylation) table and run Fisher's exact
test. Call ASM when p < alpha (thresholded at PLOT time, not here — we store the
raw p-value).

    table = [[a, b],      a/b = meth/unmeth calls among H1 reads
             [c, d]]       c/d = meth/unmeth calls among H2 reads

Naive Fisher CANNOT use noH reads (no allele label) — they are excluded. This is
on purpose: it shows how much extra information the EM/BIC method exploits (it
marginalises noH reads over alleles and even works fully unphased). Hence under
the `all_noh` treatment Fisher has no table at all -> no ASM call (recall 0).

This pipeline is independent of run_sweep.py: it imports run_sweep ONLY to reuse
its generation, treatment, seeding, and grid helpers, so it scores the SAME reads
(same seed_base) as the BIC run — a paired comparison — without modifying it.
Output: one TSV shard with columns
  truth_class, depth, noh_frac, rep, approach, fisher_n_reads, fisher_p
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# Reuse the BIC harness's generation/treatment/seed helpers (read-only import).
_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.join(_HERE, os.pardir)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

import run_sweep as rs  # noqa: E402

# MUST match run_sweep.py's --seed-base default so the generated reads are
# byte-identical to the BIC run (paired comparison). If you change one, change both.
SEED_BASE = 20260603

OUTPUT_COLUMNS = (
    "truth_class", "depth", "noh_frac", "rep", "approach",
    "fisher_n_reads", "fisher_p",
)


def _fisher_score(X, tags):
    """Naive pooled Fisher exact test on the allele x methylation 2x2 table.

    Returns (n_reads_used, p). noH reads are excluded. p is NaN (no ASM call)
    when there are no H1/H2 reads (e.g. the all_noh treatment) or when either
    allele contributes zero observed calls.
    """
    sel = (tags == "H1") | (tags == "H2")
    n_used = int(sel.sum())
    if n_used == 0:
        return 0, float("nan")

    Xs = X[sel]
    ts = tags[sel]
    h1 = ts == "H1"
    h2 = ts == "H2"
    meth = Xs == 1      # NaN compares False in both masks
    unmeth = Xs == 0

    a = int(meth[h1].sum())
    b = int(unmeth[h1].sum())
    c = int(meth[h2].sum())
    d = int(unmeth[h2].sum())

    if (a + b) == 0 or (c + d) == 0:
        return n_used, float("nan")

    p = float(fisher_exact([[a, b], [c, d]])[1])
    return n_used, p


def run_cell(region_class, depth, noh_frac, noh_idx, reps, seed_base):
    """Yield Fisher rows for all reps x approaches of one (class, depth, noH) cell.

    Mirrors run_sweep.run_cell's generation/treatment sequence exactly (same
    seed -> same reads), swapping the EM/BIC scorer for the Fisher scorer.
    """
    rows = []
    class_id = rs._CLASS_ID[region_class]
    for rep in range(reps):
        seed = rs._stable_seed(seed_base, class_id, depth, noh_idx, rep)
        rng = np.random.RandomState(seed)
        X, tags = rs._build_region(region_class, depth, noh_frac, rng)
        for approach in rs.APPROACHES:
            Xv, tv = rs._treatment(approach, X, tags, rng)
            n_used, p = _fisher_score(Xv, tv)
            rows.append({
                "truth_class": region_class,
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
    parser.add_argument("--grid", choices=rs.GRIDS, required=True)
    parser.add_argument("--task-id", type=int, required=True,
                        help="SGE array task id (1-based); selects one depth.")
    parser.add_argument("--out", type=str, required=True, help="Output shard TSV path.")
    parser.add_argument("--reps", type=int, default=None,
                        help="Regions per class per cell (default: grid-specific).")
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
        f"[fisher] grid={args.grid} task={args.task_id}/{len(depths)} depth={depth} "
        f"noH={len(nohs)} reps/class={reps} seed_base={args.seed_base}"
    )

    all_rows = []
    for noh_idx, noh_frac in rs.noh_cells(args.grid):
        for region_class in (rs.gs.SHARED, rs.gs.STATIC_ASM):
            all_rows.extend(run_cell(
                region_class, depth, noh_frac, noh_idx, reps, args.seed_base,
            ))

    out = pd.DataFrame(all_rows, columns=list(OUTPUT_COLUMNS))
    out = out.sort_values(["truth_class", "noh_frac", "rep", "approach"]).reset_index(drop=True)

    tmp = args.out + ".tmp"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, args.out)
    print(f"[fisher] Wrote {len(out)} rows to {args.out}")


if __name__ == "__main__":
    main()
