# ABOUTME: Driver for static ASM detection — Shared (M0) vs ASM (M1) BIC per window.
# ABOUTME: Single-sample input; emits a TSV of per-window BICs. No time dimension.
import argparse
import glob
import os
import warnings
from multiprocessing import Pool

import EMfunctions_ASM as em


def main():
    warnings.filterwarnings("ignore", category=FutureWarning, module="numpy")

    parser = argparse.ArgumentParser(
        description="Static ASM detection using EM + BIC on a single sample. "
                    "Computes Shared (M0) and ASM (M1) BICs per CpG window; "
                    "BICresult = 1 iff ASM beats Shared."
    )
    parser.add_argument(
        'input_path',
        help='Directory of per-region TSVs produced by '
             'bucket_reads_into_regions.py (read_format_split/<chr>/).'
    )
    parser.add_argument('--output', '-o', type=str, required=True, help='Output TSV path')
    parser.add_argument('--windowsize', type=int, default=10, help='CpG sites per window')
    parser.add_argument(
        '--window-step', type=int, default=em.DEFAULT_WINDOW_STEP,
        help=(
            'CpGs advanced between consecutive evaluated windows (default '
            f'{em.DEFAULT_WINDOW_STEP}). 1 = evaluate every in-region CpG '
            '(fully overlapping); larger = coarser, faster scan.'
        ),
    )
    parser.add_argument('--workers', type=int, default=10,
                        help='Parallel worker processes (default 10).')
    parser.add_argument('--glob-pattern', type=str, default='methylationfraction_*_.tsv',
                        help='Glob pattern for per-region TSVs.')
    parser.add_argument(
        '--asm-restarts', type=int, default=em.DEFAULT_RESTARTS,
        help=(
            'EM restarts for the ASM (M1) BIC, keep min BIC. '
            'Default 5 matches the 2026-05-18 restart calibration. Use 1 '
            'for A/B against the pre-restart baseline.'
        ),
    )
    parser.add_argument(
        '--em-tol', type=float, default=em.DEFAULT_EM_TOL,
        help=(
            'Absolute convergence tolerance for the ASM (M1) EM: stops when '
            f'abs(ll - previous_ll) < tol. Default {em.DEFAULT_EM_TOL:g}. '
            'Single fixed absolute value (see EM-tol calibration).'
        ),
    )
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_path, args.glob_pattern)))
    print(f"{len(files)} region files in {args.input_path}")
    if not files:
        raise SystemExit(f"No region files matched {args.glob_pattern} in {args.input_path}")

    print(f"EM restarts: ASM={args.asm_restarts}")
    print(f"window step: {args.window_step}")
    print(f"EM tol: {args.em_tol:g}")

    tasks = [(p, args.windowsize, args.asm_restarts, args.window_step, args.em_tol) for p in files]

    tmp_output = args.output + ".tmp"
    header_written = False
    with Pool(processes=args.workers) as pool, open(tmp_output, 'w') as f:
        for df in pool.imap_unordered(_process_region, tasks):
            df.to_csv(f, sep='\t', index=False, header=not header_written)
            header_written = True
            f.flush()
    os.replace(tmp_output, args.output)
    print(f"Output written to {args.output}")


def _process_region(task):
    path, windowsize, asm_restarts, window_step, em_tol = task
    return em.EMBIC_bin_path_ASM(
        path,
        windowsize=windowsize,
        asm_restarts=asm_restarts,
        window_step=window_step,
        tol=em_tol,
    )


if __name__ == "__main__":
    main()
