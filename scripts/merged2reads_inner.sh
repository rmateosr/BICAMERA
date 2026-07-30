#!/bin/bash
# ABOUTME: SGE inner — per-read aggregation of the merged coordmod TSV into one row per read.
# ABOUTME: Single awk pass; replaces the previous Rscript line-by-line implementation.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=2G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset
set -o pipefail

OUTPUT_FOLDER=$1
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
CHR=${chromosomes[$((SGE_TASK_ID-1))]}
echo "$CHR"

IN="${OUTPUT_FOLDER}/modkit_referenced_splitbam_mergedbychr/modkit_extract_${CHR}_merged_filtered_coordmod.tsv"
OUT_DIR="${OUTPUT_FOLDER}/read_format"
OUT="${OUT_DIR}/${CHR}_reads.tsv"
TMP="${OUT}.tmp"

mkdir -p "${OUT_DIR}"

# Input (no header, tab-separated, produced by coordfixer.py):
#   col 1  read_id
#   col 2  CpG position
#   col 3  chrom
#   col 4  modkit call_code (string: "m" methylated, "-" canonical)
#   col 5  haplotype (H1 | H2 | noH)
#
# Assumed sort: rows for the same read_id appear contiguously. coordfixer.py
# concatenates H1/H2/noH blocks; within each block modkit emits per-read CpGs
# contiguously, and a single physical read can only belong to one haplotype,
# so the contiguity assumption holds without an explicit sort.
#
# Output (header + one row per read, tab-separated):
#   readlabel  haplotype  chrom  startcoord  status
# where startcoord and status are comma-joined values from the read's CpGs.

awk -F'\t' -v OFS='\t' '
BEGIN {
    print "readlabel", "haplotype", "chrom", "startcoord", "status"
}
{
    rid = $1
    cpg = $2
    chrom = $3
    meth = ($4 == "m") ? 1 : 0
    hap = $5

    if (rid != prev_rid && prev_rid != "") {
        print prev_rid, prev_hap, prev_chrom, coords, statuses
        coords = ""
        statuses = ""
    }
    if (coords == "") {
        coords = cpg
        statuses = meth
    } else {
        coords = coords "," cpg
        statuses = statuses "," meth
    }
    prev_rid = rid
    prev_hap = hap
    prev_chrom = chrom
}
END {
    if (prev_rid != "") {
        print prev_rid, prev_hap, prev_chrom, coords, statuses
    }
}
' "${IN}" > "${TMP}"
mv "${TMP}" "${OUT}"
