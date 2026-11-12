#!/bin/bash
# ABOUTME: Cuts a fine-grid sweep TSV (bic/fisher/asms) down to the manuscript's paper grid, no re-run.
# ABOUTME: Usage: bash benchmarks/subset_paper_grid.sh <fine.tsv> <out.tsv>   (keeps depth 10..100 x noH 0.0..1.0 by 0.1)
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fine.tsv> <out.tsv>" >&2
    exit 2
fi
IN=$1
OUT=$2
[ -f "${IN}" ] || { echo "ERROR: file not found: ${IN}" >&2; exit 1; }
[ "$(head -1 "${IN}" | cut -f2,3)" = "$(printf 'depth\tnoh_frac')" ] || {
    echo "ERROR: ${IN}: expected columns 2-3 to be depth, noh_frac" >&2; exit 1; }

# Keep depth multiples of 10 (>= 10) and noH multiples of 0.10. The *100+0.5 round
# guards against a value such as 0.3 printed as 0.30000000000000004.
awk -F'\t' 'NR==1{print; next}
    {d=$2+0; m=int($3*100+0.5); if (d>=10 && d%10==0 && m%10==0) print}' "${IN}" > "${OUT}.tmp"
mv "${OUT}.tmp" "${OUT}"

echo "Wrote $(($(wc -l < "${OUT}") - 1)) rows to ${OUT}"
