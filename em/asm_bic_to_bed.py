#!/usr/bin/env python3
# ABOUTME: Convert per-window ASM BIC TSVs (<chr>_ASM.tsv) into whole-genome BEDs of ASM windows for IGV.
# ABOUTME: One run from -o PREFIX writes <PREFIX>_scored.bed, <PREFIX>_merged.bed, <PREFIX>_all_colored.bed.

# Each row of <chr>_ASM.tsv is one EM window. `coord_evaluated` is the 0-based
# ref_position of the first CpG in the window; `window_last_coord` is the 0-based
# ref_position of the last CpG. A window is ASM iff `BICresult == 1` (ASM/M1 beat
# Shared/M0). The ASM margin = BICmiddlecomp - BICsinglecomp is < 0 exactly for
# ASM windows (lower BIC wins, matching BICresult); sentinel low-cov windows -> BICresult==0.
#
# BED is 0-based half-open, so to make the interval cover the last CpG base we emit
# end = window_last_coord + 1. Three outputs per run:
#   <PREFIX>_scored.bed      : one interval per ASM window, BED4 with a 4th column
#                              = ASM margin (middle-single; more negative = stronger ASM). ASM only.
#   <PREFIX>_merged.bed      : overlapping/bookending ASM windows (consecutive
#                              windows step by < windowsize CpGs inside a region)
#                              collapsed into contiguous BED3 blocks — cleaner in
#                              IGV than stacked boxes. ASM only.
#   <PREFIX>_all_colored.bed : EVERY analyzed window (ASM and non-ASM), not merged,
#                              BED9 with itemRgb so IGV colors ASM red and non-ASM
#                              blue. col4 name = signed margin (middle-single, <0=ASM);
#                              col5 score = |margin| clamped to BED's 0-1000 int (BED
#                              score must be non-negative; the ASM call is carried by
#                              color). Carries a `track ... itemRgb="On"` header line.

import argparse
import glob
import os
import sys


def chrom_sort_key(chrom):
    """Natural order: chr1..chr22, chrX, chrY, chrM, then anything else."""
    name = chrom[3:] if chrom.startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), name)
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name in order:
        return (0, order[name], name)
    return (1, 0, chrom)


def chrom_from_filename(path):
    """chr15_ASM.tsv -> chr15."""
    base = os.path.basename(path)
    return base[: -len("_ASM.tsv")] if base.endswith("_ASM.tsv") else base.split(".")[0]


def parse_asm_tsv(path, chrom):
    """Yield (chrom, start, end, score, is_asm) for EVERY analyzed window in one TSV.

    score is the ASM margin = BICmiddlecomp - BICsinglecomp (more negative = stronger
    ASM; < 0 exactly when the window is ASM). is_asm is True iff BICresult == 1.
    """
    intervals = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            i_first = header.index("coord_evaluated")
            i_last = header.index("window_last_coord")
            i_res = header.index("BICresult")
            i_m0 = header.index("BICsinglecomp")
            i_m1 = header.index("BICmiddlecomp")
        except ValueError as e:
            sys.exit(f"ERROR: {path} missing expected column: {e}")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            start = int(fields[i_first])
            end = int(fields[i_last]) + 1  # half-open, include last CpG base
            score = float(fields[i_m1]) - float(fields[i_m0])  # BICmiddlecomp - BICsinglecomp; <0 = ASM
            is_asm = fields[i_res] == "1"
            intervals.append((chrom, start, end, score, is_asm))
    return intervals


def merge_intervals(intervals):
    """Collapse overlapping/bookending intervals per chrom. Drops scores."""
    by_chrom = {}
    for chrom, start, end, *_ in intervals:
        by_chrom.setdefault(chrom, []).append((start, end))
    merged = []
    for chrom in by_chrom:
        spans = sorted(by_chrom[chrom])
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e:  # overlap or bookend (half-open touch)
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e, None))
                cur_s, cur_e = s, e
        merged.append((chrom, cur_s, cur_e, None))
    return merged


# itemRgb colors for the all-windows track (RColorBrewer Set1).
ASM_RGB = "228,26,28"      # red  (#E41A1C) — ASM (BICresult==1)
NON_ASM_RGB = "55,126,184"  # blue (#377EB8) — non-ASM (BICresult==0)


def main():
    parser = argparse.ArgumentParser(
        description="Whole-genome BEDs of ASM windows from per-chr <chr>_ASM.tsv."
    )
    parser.add_argument(
        "asm_dir",
        help="Directory containing <chr>_ASM.tsv files (e.g. <out>/ASM).",
    )
    parser.add_argument(
        "-o", "--out-prefix", required=True,
        help="Output prefix; writes <PREFIX>_scored.bed (ASM-only, BED4), "
             "<PREFIX>_merged.bed (ASM-only contiguous blocks, BED3), and "
             "<PREFIX>_all_colored.bed (every window, BED9, ASM red / non-ASM blue).",
    )
    parser.add_argument(
        "--glob", default="*_ASM.tsv",
        help="Filename pattern inside asm_dir (default '*_ASM.tsv').",
    )
    args = parser.parse_args()

    tsvs = sorted(glob.glob(os.path.join(args.asm_dir, args.glob)))
    if not tsvs:
        sys.exit(f"ERROR: no files matching {args.glob} under {args.asm_dir}")

    all_windows = []   # every analyzed window: (chrom, start, end, score, is_asm)
    for tsv in tsvs:
        chrom = chrom_from_filename(tsv)
        ivs = parse_asm_tsv(tsv, chrom)
        all_windows.extend(ivs)
        n_asm = sum(1 for x in ivs if x[4])
        print(f"{chrom}: {len(ivs)} windows ({n_asm} ASM)", file=sys.stderr)

    win_sort = lambda x: (chrom_sort_key(x[0]), x[1], x[2])
    asm_windows = [x for x in all_windows if x[4]]

    # 1) Per-window scored BED, ASM only (BED4).
    scored = sorted(asm_windows, key=win_sort)
    scored_path = f"{args.out_prefix}_scored.bed"
    with open(scored_path + ".tmp", "w") as f:
        for chrom, start, end, score, _ in scored:
            f.write(f"{chrom}\t{start}\t{end}\t{score:.2f}\n")
    os.replace(scored_path + ".tmp", scored_path)
    print(f"wrote {len(scored)} scored ASM windows to {scored_path}", file=sys.stderr)

    # 2) Merged contiguous-block BED, ASM only (BED3).
    merged = merge_intervals(asm_windows)
    merged.sort(key=win_sort)
    merged_path = f"{args.out_prefix}_merged.bed"
    with open(merged_path + ".tmp", "w") as f:
        for chrom, start, end, _ in merged:
            f.write(f"{chrom}\t{start}\t{end}\n")
    os.replace(merged_path + ".tmp", merged_path)
    print(f"wrote {len(merged)} merged ASM blocks to {merged_path}", file=sys.stderr)

    # 3) All analyzed windows, colored by ASM/non-ASM (BED9 + itemRgb).
    every = sorted(all_windows, key=win_sort)
    all_path = f"{args.out_prefix}_all_colored.bed"
    track = (
        f'track name="{os.path.basename(args.out_prefix)}_all" '
        f'description="per-window ASM (red) vs non-ASM (blue); name/score = BIC margin" '
        f'itemRgb="On"\n'
    )
    with open(all_path + ".tmp", "w") as f:
        f.write(track)
        for chrom, start, end, score, is_asm in every:
            rgb = ASM_RGB if is_asm else NON_ASM_RGB
            # col5 score must be a non-negative int in [0, 1000]; use |margin| so display
            # intensity tracks evidence strength (col4 holds the signed margin; the
            # ASM/non-ASM call is carried by the RGB color).
            score_int = max(0, min(1000, round(abs(score))))
            f.write(
                f"{chrom}\t{start}\t{end}\t{score:.2f}\t{score_int}\t.\t"
                f"{start}\t{end}\t{rgb}\n"
            )
    os.replace(all_path + ".tmp", all_path)
    n_asm = sum(1 for x in every if x[4])
    print(
        f"wrote {len(every)} windows ({n_asm} ASM, {len(every) - n_asm} non-ASM) "
        f"to {all_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
