#!/bin/bash
# ABOUTME: One array task of a benchmark sweep: runs the sweep script for task SGE_TASK_ID into one shard.
# ABOUTME: Queued by run_benchmarks.sh. Args: <script under benchmarks/> <shard_dir> <label> [script args].
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
set -euo pipefail

# Repo root: exported by run_benchmarks.sh; SGE_O_WORKDIR / $(pwd) cover a task run by hand.
BICAMERA_ROOT="${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}"
. "${BICAMERA_ROOT}/config/site.sh"
bicamera_require python

# One BLAS thread: the matrices are tiny, and an unpinned OpenBLAS starts one
# thread per core, which can exceed the job's process or memory limits.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

SCRIPT=$1
SHARD_DIR=$2
LABEL=$3
shift 3
: "${SGE_TASK_ID:?SGE_TASK_ID not set — run through run_benchmarks.sh or set it by hand}"

mkdir -p "${SHARD_DIR}"
echo "=== ${LABEL} task ${SGE_TASK_ID}: ${SCRIPT} $* ==="
"${PYTHON}" "${BICAMERA_ROOT}/benchmarks/${SCRIPT}" \
    --task-id "${SGE_TASK_ID}" \
    --out "${SHARD_DIR}/${LABEL}_task_${SGE_TASK_ID}.tsv" \
    "$@"
echo "=== ${LABEL} task ${SGE_TASK_ID} done ==="
