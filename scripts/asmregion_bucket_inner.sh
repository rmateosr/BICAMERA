#!/bin/bash
# ABOUTME: SGE inner — bucket a chromosome's read_format rows into per-region TSVs.
# ABOUTME: Replaces the divide step; boundaries come from potential_asm_regions/<chr>.bed.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=4G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load python/3.12.0

OUTPUT_FOLDER=$1
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
CHR=${chromosomes[$((SGE_TASK_ID-1))]}
echo "$CHR"

python3.12 scripts/bucket_reads_into_regions.py "${CHR}" "${OUTPUT_FOLDER}"
