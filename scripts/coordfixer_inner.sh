#!/bin/bash
# ABOUTME: Per-chr coordinate fix + haplotype join over the filtered modkit extract.
# ABOUTME: Thin SGE wrapper around scripts/coordfixer.py.
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=2G
set -xv
set -o errexit
set -o nounset

# Repo root: exported by run_pipeline.sh. SGE_O_WORKDIR / $(pwd) cover a stage
# invoked by hand. site.sh resolves the tools; defaults.sh supplies the tunables.
# This stage previously called `python3.12` with no module load at all, working
# only because qsub -V happened to carry a PATH containing that exact binary name.
BICAMERA_ROOT="${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}"
. "${BICAMERA_ROOT}/config/site.sh"
. "${BICAMERA_ROOT}/config/defaults.sh"
bicamera_require python

# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
CHR=${chromosomes[$((SGE_TASK_ID-1))]}
OUTPUT_FOLDER=$1

"${PYTHON}" "${BICAMERA_ROOT}/scripts/coordfixer.py" "${CHR}" "${OUTPUT_FOLDER}"



