#!/bin/bash
# ABOUTME: Prep an ONT modBAM for the 5mC ASM pipeline: strip 5hmC (C+h) then keep MAPQ>=THRESH.
# ABOUTME: Reusable across ONT samples. Idempotent, atomic; prints the final BAM path on stdout.
set -euo pipefail

usage() {
  cat >&2 <<EOF
Usage: $0 [-q MAPQ] [-t THREADS] [-m MODKIT_BIN] <in.bam>

  <in.bam>      ONT modBAM with C+m and C+h MM/ML tags (e.g. dorado sup).
  -q MAPQ       minimum MAPQ to keep (default: 10; samtools -q is inclusive, so 10 == MAPQ>=10).
  -t THREADS    threads for modkit/samtools (default: 1; does not change output, only speed).
  -m MODKIT_BIN modkit binary (default: resolved from PATH by config/site.sh).

Two steps, both skipped if their output already exists:
  1. modkit adjust-mods --ignore h  ->  <in>.5mc.bam        (5hmC removed; 5mC-vs-canonical only)
  2. samtools view -q MAPQ          ->  <in>.5mc.mapq<Q>.bam (the pipeline input)

Why strip 5hmC: otherwise 'modkit extract calls' emits 'h' calls that the pipeline
filter silently drops AND the canonical-vs-5mC threshold is estimated 3-way — making
ONT non-comparable to 5mC-only HiFi data. See docs/input_requirements.md.

Prints the final BAM path to stdout (so callers can capture it).
EOF
  exit 1
}

export BICAMERA_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Resolves MODKIT_BIN and SAMTOOLS from PATH, applying a site profile if one
# matches. -m still overrides modkit explicitly.
. "${BICAMERA_ROOT}/config/site.sh"

MAPQ=10
THREADS=1
while getopts "q:t:m:h" opt; do
  case "$opt" in
    q) MAPQ=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    m) MODKIT_BIN=$OPTARG ;;
    h|*) usage ;;
  esac
done
shift $((OPTIND - 1))
[[ $# -eq 1 ]] || usage
IN=$1

[[ -f $IN ]] || { echo "ERROR: input BAM not found: $IN" >&2; exit 1; }
bicamera_require modkit samtools

base="${IN%.bam}"
NOH="${base}.5mc.bam"
OUT="${base}.5mc.mapq${MAPQ}.bam"

if [[ -f $OUT && -f ${OUT}.bai ]]; then
  echo "SKIP prep (exists): $OUT" >&2
  echo "$OUT"; exit 0
fi

if [[ ! -f $NOH || ! -f ${NOH}.bai ]]; then
  echo "[$(date +%T)] strip 5hmC: adjust-mods --ignore h -> $NOH" >&2
  "$MODKIT_BIN" adjust-mods --ignore h -t "$THREADS" "$IN" "${NOH}.tmp"
  mv "${NOH}.tmp" "$NOH"
  "$SAMTOOLS" index -@ "$THREADS" "$NOH"
fi

echo "[$(date +%T)] keep MAPQ>=${MAPQ} -> $OUT" >&2
"$SAMTOOLS" view -@ "$THREADS" -b -q "$MAPQ" "$NOH" -o "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"
"$SAMTOOLS" index -@ "$THREADS" "$OUT"
echo "[$(date +%T)] done. reads: $("$SAMTOOLS" view -c "$OUT")" >&2

echo "$OUT"
