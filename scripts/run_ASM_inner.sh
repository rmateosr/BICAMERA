#!/bin/bash
# ABOUTME: Per-chromosome SGE job — runs Shared+ASM EM+BIC over per-region TSVs for one sample.
# ABOUTME: Invokes em/BIC_algorithm_ASM.py with 10 workers; no time dimension, no multi-sample pairing.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=2G
#$ -pe def_slot 10
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load python/3.12.0

# Pin each worker's BLAS backend to 1 thread.
# Rationale: with WORKERS=10 and BLAS defaults (one pool per worker, each ~= physical cores),
# the process would spawn >>10 threads and contend with itself.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Args: OUT_DIR SAMPLE_OUT  (chromosome derived from SGE_TASK_ID)
OUT_DIR=$1
SAMPLE_OUT=$2
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
CHR=${chromosomes[$((SGE_TASK_ID-1))]}

WINDOWSIZE=${WINDOWSIZE:-10}
WORKERS=${WORKERS:-10}
WINDOW_STEP=${WINDOW_STEP:-5}

RESTART_FLAGS=()
if [ -n "${ASM_RESTARTS:-}" ]; then
    RESTART_FLAGS+=(--asm-restarts "${ASM_RESTARTS}")
fi

EM_TOL_FLAGS=()
if [ -n "${EM_TOL:-}" ]; then
    EM_TOL_FLAGS+=(--em-tol "${EM_TOL}")
fi

mkdir -p "${OUT_DIR}"

REGION_DIR="${SAMPLE_OUT}/read_format_split/${CHR}"
if [ ! -d "${REGION_DIR}" ]; then
    echo "ERROR: missing per-chr region dir: ${REGION_DIR}" >&2
    exit 1
fi

# Potential-ASM-region route: a chromosome with no intermediate-methylation
# regions yields an empty region dir. That is a valid result, not an error —
# skip the EM rather than failing the job.
shopt -s nullglob
region_files=("${REGION_DIR}"/methylationfraction_*_.tsv)
shopt -u nullglob
if [ "${#region_files[@]}" -eq 0 ]; then
    echo "NOTE: no potential-ASM region files in ${REGION_DIR}; skipping EM for ${CHR}." >&2
    exit 0
fi

# ASM_DRIVER is exported by the runner; falls back to the in-repo path for direct invocations.
ASM_DRIVER="${ASM_DRIVER:-${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}/em/BIC_algorithm_ASM.py}"

python3 "${ASM_DRIVER}" \
    "${REGION_DIR}" \
    --windowsize "${WINDOWSIZE}" \
    --window-step "${WINDOW_STEP}" \
    --workers "${WORKERS}" \
    "${RESTART_FLAGS[@]}" \
    "${EM_TOL_FLAGS[@]}" \
    --output "${OUT_DIR}/${CHR}_ASM.tsv"
