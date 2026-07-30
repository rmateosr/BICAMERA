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
  -m MODKIT_BIN modkit binary (default: <script_dir>/bin/modkit).

Two steps, both skipped if their output already exists:
  1. modkit adjust-mods --ignore h  ->  <in>.5mc.bam        (5hmC removed; 5mC-vs-canonical only)
  2. samtools view -q MAPQ          ->  <in>.5mc.mapq<Q>.bam (the pipeline input)

Why strip 5hmC: otherwise 'modkit extract calls' emits 'h' calls that the pipeline
filter silently drops AND the canonical-vs-5mC threshold is estimated 3-way — making
ONT non-comparable to 5mC-only HiFi data. See HG002_ONT_giab2023.05/PROVENANCE.md.

Prints the final BAM path to stdout (so callers can capture it).
EOF
  exit 1
}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MAPQ=10
THREADS=1
MODKIT_BIN="${SCRIPT_DIR}/bin/modkit"
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

[[ -f $IN ]]           || { echo "ERROR: input BAM not found: $IN" >&2; exit 1; }
[[ -x $MODKIT_BIN ]]   || { echo "ERROR: modkit not executable: $MODKIT_BIN" >&2; exit 1; }

module use /usr/local/package/modulefiles/ >/dev/null 2>&1 || true
module load samtools/1.19 >/dev/null 2>&1 || true

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
  samtools index -@ "$THREADS" "$NOH"
fi

echo "[$(date +%T)] keep MAPQ>=${MAPQ} -> $OUT" >&2
samtools view -@ "$THREADS" -b -q "$MAPQ" "$NOH" -o "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"
samtools index -@ "$THREADS" "$OUT"
echo "[$(date +%T)] done. reads: $(samtools view -c "$OUT")" >&2

echo "$OUT"
