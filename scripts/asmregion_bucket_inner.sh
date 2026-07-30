#!/bin/bash
# ABOUTME: SGE inner — bucket a chromosome's read_format rows into per-region TSVs.
# ABOUTME: Replaces the divide step; boundaries come from potential_asm_regions/<chr>.bed.
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=4G
#$ -pe def_slot 1
set -xv
set -o errexit
set -o nounset

# Repo root: exported by run_pipeline.sh. SGE_O_WORKDIR / $(pwd) cover a stage
# invoked by hand. site.sh resolves the tools; defaults.sh supplies the tunables.
BICAMERA_ROOT="${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}"
. "${BICAMERA_ROOT}/config/site.sh"
. "${BICAMERA_ROOT}/config/defaults.sh"
bicamera_require python

OUTPUT_FOLDER=$1
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
CHR=${chromosomes[$((SGE_TASK_ID-1))]}

# MIN_REGION_CPGS is passed explicitly rather than read from the environment
# inside the Python, so the threshold cannot silently differ from the launcher's.
"${PYTHON}" "${BICAMERA_ROOT}/scripts/bucket_reads_into_regions.py" \
    "${CHR}" "${OUTPUT_FOLDER}" --min-region-cpgs "${MIN_REGION_CPGS}"
