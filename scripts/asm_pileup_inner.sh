#!/bin/bash
# ABOUTME: SGE inner — per-chr modkit pileup of the whole BAM into per-CpG bedMethyl.
# ABOUTME: Population-average methylation; feeds the potential-ASM-region selector.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=8G
#$ -pe def_slot 4
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load samtools/1.19
module load apptainer/

OUTPUT_FOLDER=$1
REF=$2
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
chr=${chromosomes[$((SGE_TASK_ID-1))]}

# Same modkit-binary resolution rationale as modkitinner_withreference_bychrandhap.sh:
# under SGE $0 is the spool path, so resolve via $BICAMERA_ROOT / $SGE_O_WORKDIR.
MODKIT_BIN=${MODKIT_BIN:-${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}/bin/modkit}

mkdir -p "${OUTPUT_FOLDER}/pileup"
OUT="${OUTPUT_FOLDER}/pileup/${chr}.bedmethyl"
TMP="${OUT}.tmp"

# Population-average per-CpG methylation across ALL reads (H1/H2/noH together).
# --cpg --combine-strands -> one row per CpG dyad; --ignore h folds any 5hmC into
# canonical (harmless: prep already stripped 5hmC). --threads 4 matches the
# `#$ -pe def_slot 4` allocation: after the splitbam stage was deleted, the
# pileup->regions discovery leg is the longer of the two legs feeding modkitref,
# so pileup gates modkitref's start and is worth threading. Output is
# thread-count-independent (verified byte-identical to --threads 1 on chr15).
# bedMethyl: col 10 = Nvalid_cov, col 11 = fraction_modified (percent).
${MODKIT_BIN} pileup --threads 4 --cpg --combine-strands --ignore h \
    --ref "${REF}" \
    --region "${chr}" \
    "${OUTPUT_FOLDER}/phased_output.bam" \
    "${TMP}"
mv "${TMP}" "${OUT}"
