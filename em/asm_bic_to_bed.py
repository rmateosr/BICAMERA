#!/usr/bin/env python3
# ABOUTME: Convert per-window ASM BIC TSVs (<chr>_ASM.tsv) into whole-genome BEDs of ASM windows for IGV.
# ABOUTME: One run from -o PREFIX writes <PREFIX>_scored.bed, <PREFIX>_merged.bed, <PREFIX>_all_colored.bed.

# Each row of <chr>_ASM.tsv is one EM window. `coord_evaluated` is the 0-based
# ref_position of the first CpG in the window; `window_last_coord` is the 0-based
# ref_position of the last CpG. A window is ASM iff `BICresult == 1` (ASM/M1 beat
# Shared/M0).
#
# Signed evidence is reported as delta BIC, in the same direction the paper defines:
#
#     dBIC = BICsinglecomp - BICmiddlecomp = BIC_Shared - BIC_ASM
#
# so dBIC > 0 favours ASM and larger means stronger. This agrees with BICresult by
# construction: BICresult = int(BICmiddlecomp < BICsinglecomp), which is exactly
# dBIC > 0. Exact ties give dBIC == 0 and BICresult == 0 (non-ASM), and the
# low-coverage sentinel windows likewise come out non-ASM.
#
# BED is 0-based half-open, so to make the interval cover the last CpG base we emit
# end = window_last_coord + 1. Three outputs per run:
#   <PREFIX>_scored.bed      : one interval per ASM window, BED4 with a 4th column
#                              = dBIC (> 0; larger = stronger ASM). ASM only.
#   <PREFIX>_merged.bed      : overlapping/bookending ASM windows (consecutive
#                              windows step by < windowsize CpGs inside a region)
#                              collapsed into contiguous BED3 blocks — cleaner in
#                              IGV than stacked boxes. ASM only.
#   <PREFIX>_all_colored.bed : EVERY analyzed window (ASM and non-ASM), not merged,
#                              BED9 with itemRgb so IGV colors ASM red and non-ASM
#                              blue. col4 name = signed dBIC (> 0 = ASM);
#                              col5 score = |dBIC| clamped to BED's 0-1000 int. The
#                              absolute value is REQUIRED: non-ASM windows have
#                              dBIC < 0 and BED col5 must be a non-negative int, so
#                              dropping abs() would clamp every non-ASM window to 0.
#                              The ASM call itself is carried by the color.
#                              Carries a `track ... itemRgb="On"` header line.

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
    """Yield (chrom, start, end, delta_bic, is_asm) for EVERY analyzed window in one TSV.

    delta_bic = BICsinglecomp - BICmiddlecomp = BIC_Shared - BIC_ASM, so it is > 0
    exactly when the window is ASM (larger = stronger). is_asm is True iff
    BICresult == 1, which is the same condition by construction.
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
        # The driver emits deltaBIC directly. Fall back to deriving it for TSVs
        # written before that column existed, so older runs still convert.
        i_delta = header.index("deltaBIC") if "deltaBIC" in header else None
        for line in f:
            fields = line.rstrip("\n").split("\t")
            start = int(fields[i_first])
            end = int(fields[i_last]) + 1  # half-open, include last CpG base
            if i_delta is not None:
                delta_bic = float(fields[i_delta])
            else:
                # BICsinglecomp - BICmiddlecomp = BIC_Shared - BIC_ASM; > 0 = ASM.
                delta_bic = float(fields[i_m0]) - float(fields[i_m1])
            is_asm = fields[i_res] == "1"
            intervals.append((chrom, start, end, delta_bic, is_asm))
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

# Print format for the signed dBIC written to BED col4.
#
# Significant digits, not fixed decimals, because dBIC spans several orders of
# magnitude: on chr22 of a 48x HiFi genome it runs from 0.00033 to 221.35. Any fixed
# %.Nf renders everything below its resolution as a signless "0.00...", which makes a
# red ASM window read as a tie and breaks the "dBIC > 0 for ASM" invariant in the
# emitted file (never in the arithmetic). At 2 decimals that hit 3 of chr22's 8,251
# ASM windows; picking a deeper fixed precision only moves the cliff. %g has no such
# cliff — a nonzero dBIC never prints as zero at any magnitude.
#
# The per-window TSV remains the full-precision source of truth; this is the BED's
# display/interchange label.
DBIC_FMT = ".6g"


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

    all_windows = []   # every analyzed window: (chrom, start, end, delta_bic, is_asm)
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
        for chrom, start, end, delta_bic, _ in scored:
            f.write(f"{chrom}\t{start}\t{end}\t{delta_bic:{DBIC_FMT}}\n")
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
        f'description="per-window ASM (red) vs non-ASM (blue); '
        f'name = dBIC (BIC_Shared - BIC_ASM, >0 = ASM), score = |dBIC|" '
        f'itemRgb="On"\n'
    )
    with open(all_path + ".tmp", "w") as f:
        f.write(track)
        for chrom, start, end, delta_bic, is_asm in every:
            rgb = ASM_RGB if is_asm else NON_ASM_RGB
            # col5 must be a non-negative int in [0, 1000], so take |dBIC| and let
            # display intensity track evidence strength. abs() is load-bearing, NOT
            # redundant now that dBIC > 0 means ASM: non-ASM windows have dBIC < 0,
            # and without abs() the max(0, ...) below would clamp them all to 0.
            # col4 keeps the signed dBIC; the ASM/non-ASM call is carried by the color.
            score_int = max(0, min(1000, round(abs(delta_bic))))
            f.write(
                f"{chrom}\t{start}\t{end}\t{delta_bic:{DBIC_FMT}}\t{score_int}\t.\t"
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
