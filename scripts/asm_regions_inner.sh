#!/bin/bash
# ABOUTME: SGE inner — select intermediate-methylation CpGs and merge into potential-ASM regions.
# ABOUTME: bedMethyl -> (cov>=MIN_COV & FRAC_LOW<=frac<=FRAC_HIGH) -> bedtools merge -d MERGE_DIST.
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=2G
set -xv
set -o errexit
set -o nounset
set -o pipefail

# Repo root: exported by run_pipeline.sh. SGE_O_WORKDIR / $(pwd) cover a stage
# invoked by hand. site.sh resolves the tools; defaults.sh supplies the tunables.
BICAMERA_ROOT="${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}"
. "${BICAMERA_ROOT}/config/site.sh"
. "${BICAMERA_ROOT}/config/defaults.sh"
bicamera_require bedtools

OUTPUT_FOLDER=$1
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
chr=${chromosomes[$((SGE_TASK_ID-1))]}

# MIN_COV / FRAC_LOW / FRAC_HIGH / MERGE_DIST come from config/defaults.sh, which
# the launcher sources too — one definition, so a lost environment cannot make this
# stage select CpGs at thresholds the run was not launched with.
# See docs/region_preselection.md for the FRAC band rationale.

PILEUP="${OUTPUT_FOLDER}/pileup/${chr}.bedmethyl"
OUT_DIR="${OUTPUT_FOLDER}/potential_asm_regions"
OUT="${OUT_DIR}/${chr}.bed"
TMP="${OUT}.tmp"
mkdir -p "${OUT_DIR}"

if [ ! -s "${PILEUP}" ]; then
    echo "NOTE: empty/missing pileup for ${chr} (${PILEUP}); writing empty region BED." >&2
    : > "${TMP}"
    mv "${TMP}" "${OUT}"
    exit 0
fi

# bedMethyl: col 10 = Nvalid_cov, col 11 = fraction_modified (percent).
# Select intermediate-methylation CpGs with enough coverage, then merge CpGs
# within MERGE_DIST bp into one potential-ASM region. awk emits the chrom/start/end
# (cols 1-3); sort is required by bedtools merge.
#
# +1 on the merged region END only (the final awk). pileup runs --combine-strands, so a
# CpG dyad is reported at its C and bedMethyl end = C+1; the merged region end therefore
# lands on the terminal CpG's G (= C+1). Downstream `modkit extract calls` runs per-strand
# (NOT --combine-strands) and reports that terminal CpG's -strand 5mC at the G = the
# region's half-open end, so --include-bed EXCLUDES it and ~half the terminal CpG's reads
# (the -strand ones) are dropped. coordfixer's -strand pos-=1 cannot recover them: it runs
# after extraction, only on rows that already survived --include-bed. Extending the end by
# exactly one base pulls the terminal G inside --include-bed (coordfixer then maps it back
# onto the C). It must be +1 and no more: the next CpG's C sits at end+1 in the CGCG
# adjacency case, which the half-open --include-bed end then EXCLUDES; +2 would admit it
# and leak a methylation call from outside the region.
awk -F'\t' -v OFS='\t' -v cov="${MIN_COV}" -v lo="${FRAC_LOW}" -v hi="${FRAC_HIGH}" '
    $10 >= cov && $11 >= lo*100 && $11 <= hi*100 { print $1, $2, $3 }
' "${PILEUP}" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | "${BEDTOOLS}" merge -d "${MERGE_DIST}" \
  | awk -F'\t' -v OFS='\t' '{ print $1, $2, $3 + 1 }' \
  > "${TMP}"
mv "${TMP}" "${OUT}"

echo "${chr}: $(wc -l < "${OUT}") potential-ASM regions"
