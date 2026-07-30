#!/bin/bash
# ABOUTME: Per-chr filter of the raw modkit extract down to valid m/- CpG calls.
# ABOUTME: Pure awk on the TSV; writes atomically then deletes the raw extract.
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=2G
set -xv
set -o errexit
set -o nounset

# No config/site.sh here on purpose: this stage runs awk, mv and rm only. It
# previously loaded samtools and apptainer and invoked neither, which under
# `set -o errexit` made a renamed module kill a job that never needed the tool.
OUTPUT_FOLDER=$1
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
chr=${chromosomes[$((SGE_TASK_ID-1))]}
# Output dir was created by the (now-removed) outer; create it here so the task is self-contained.
mkdir -p "${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap_filtered"
# Keep every valid m/- call (fixed 0.5 binarization via call_code). The `fail`
# column ($18) is intentionally NOT consulted: modkit runs with --no-filtering, so
# all calls pass and we drop nothing on confidence (see modkitinner comment).
#
# De-dup: the raw extract has exactly one consumer (this stage) and is a near-byte
# duplicate of the filtered output (the predicate drops almost nothing — chr1 WG ONT:
# 263 of 12.7 GB). Write the filtered file atomically (.tmp -> mv), then delete the
# raw so only one copy persists (halves the largest TSV stage's persistent footprint).
# Re-run cost: regenerating the filtered file alone now requires re-running modkitref
# (the raw is gone). Peak transient disk is unchanged — filter holds on the whole
# modkitref array, so all raw extracts still coexist at the array boundary.
RAW="${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap/modkit_extract_${chr}.tsv"
OUT="${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap_filtered/modkit_extract_${chr}_filtered.tsv"
awk -F'\t' '$3 != "-1" && ($12 == "m" || $12 == "-") {print}' "${RAW}" > "${OUT}.tmp"
mv "${OUT}.tmp" "${OUT}"
rm -f "${RAW}"



