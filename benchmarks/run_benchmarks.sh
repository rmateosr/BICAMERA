#!/bin/bash
# ABOUTME: Runs one simulation benchmark sweep as an array of tasks plus a merge, under SGE or locally.
# ABOUTME: Usage: bash benchmarks/run_benchmarks.sh <sweep> <OUT_DIR> [extra args for the sweep script]
set -euo pipefail

# ============================================================
# Usage
# ============================================================
# bash benchmarks/run_benchmarks.sh <sweep> <OUT_DIR> [args ...]
#
# <sweep>   one of
#   bic | fisher | asms          depth x untagged-fraction grid, one task per depth.
#                                GRID=paper (default) | fine | smoke.
#   ncpg_bic | ncpg_fisher | ncpg_asms
#                                partial-ASM sweep, one task per n_asm = 0..10.
#   restarts                     EM restart calibration, one task per depth 30/50/100.
# <OUT_DIR> created if missing. Shards go to <OUT_DIR>/shards, logs to <OUT_DIR>/log,
#           the merged table to <OUT_DIR>/<label>.tsv (label = sweep, plus _<GRID>
#           for the grid sweeps, e.g. bic_paper.tsv).
# args      passed through to every task's Python script, e.g. --reps 10 for a quick
#           check. Rows depend only on (class or n_asm, depth, noH, rep), so --reps N
#           reproduces reps 0..N-1 of a full run exactly.
#
# Optional env:
#   GRID                paper | fine | smoke (bic/fisher/asms only; default paper)
#   BENCH_VMEM          per-task memory for the SGE backend (default 2G; tasks peak near 0.9 GB)
#   BICAMERA_SCHEDULER  sge | local (default: sge if qsub exists)
#   BICAMERA_SITE, PYTHON   as for run_pipeline.sh (see config/site.sh)
# ============================================================

if [ "$#" -lt 2 ]; then
    sed -n '/^# bash benchmarks/,/^# =====/p' "$0" | sed '$d;s/^# \{0,1\}//' >&2
    exit 2
fi
SWEEP=$1
OUT_DIR="$(realpath -m "$2")"
shift 2

export BICAMERA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
. "${BICAMERA_ROOT}/config/site.sh"
bicamera_require python

# Scheduler backend, chosen exactly as in run_pipeline.sh.
if [ -z "${BICAMERA_SCHEDULER+set}" ]; then
    if command -v qsub >/dev/null 2>&1; then
        BICAMERA_SCHEDULER=sge
    else
        BICAMERA_SCHEDULER=local
    fi
fi
export BICAMERA_SCHEDULER
_bicamera_backend="${BICAMERA_ROOT}/scheduler/${BICAMERA_SCHEDULER}.sh"
if [ ! -f "${_bicamera_backend}" ]; then
    echo "ERROR: unknown BICAMERA_SCHEDULER='${BICAMERA_SCHEDULER}' (no ${_bicamera_backend})" >&2
    exit 1
fi
. "${_bicamera_backend}"
unset _bicamera_backend

GRID="${GRID:-paper}"
BENCH_VMEM="${BENCH_VMEM:-2G}"

# Number of tasks for a grid sweep = number of depths in that grid, read from
# run_sweep.grid() so it cannot drift from what the tasks index.
grid_n_tasks() {
    (cd "${BICAMERA_ROOT}/benchmarks" &&
     OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     "${PYTHON}" -c 'import sys, run_sweep as rs; print(len(rs.grid(sys.argv[1])[0]))' "$1")
}

SWEEP_ARGS=()
case "${SWEEP}" in
    bic|fisher|asms)
        case "${SWEEP}" in
            bic)    SCRIPT=run_sweep.py ;;
            fisher) SCRIPT=fisher/run_fisher_sweep.py ;;
            asms)   SCRIPT=asms/run_asms_sweep.py ;;
        esac
        N=$(grid_n_tasks "${GRID}")
        SWEEP_ARGS=(--grid "${GRID}")
        LABEL="${SWEEP}_${GRID}" ;;
    ncpg_bic|ncpg_fisher|ncpg_asms)
        SCRIPT="n_asm_cpgs/run_${SWEEP}.py"
        N=11                                  # n_asm = task_id - 1 = 0..10
        LABEL="${SWEEP}" ;;
    restarts)
        SCRIPT=restarts/run_restart_sweep.py
        N=3                                   # the script's default --depths 30 50 100
        SWEEP_ARGS=(--workers 1)              # one slot per task; results do not depend on workers
        LABEL="${SWEEP}" ;;
    *)
        echo "ERROR: unknown sweep '${SWEEP}'" >&2
        echo "       expected bic, fisher, asms, ncpg_bic, ncpg_fisher, ncpg_asms or restarts" >&2
        exit 2 ;;
esac

export OUT_DIR
RUN_LABEL=bench                               # job names: bench_<label>, bench_<label>_merge
SHARD_DIR="${OUT_DIR}/shards"
MERGED="${OUT_DIR}/${LABEL}.tsv"
mkdir -p "${OUT_DIR}/log" "${SHARD_DIR}"

# Tasks run with cwd = the repo root (#$ -cwd in the task scripts).
pushd "${BICAMERA_ROOT}" >/dev/null
sweep_jid=$(bicamera_submit_array "${LABEL}" "${N}" "" \
    "s_vmem=${BENCH_VMEM},mem_req=${BENCH_VMEM}" \
    benchmarks/sweep_task.sh "${SCRIPT}" "${SHARD_DIR}" "${LABEL}" \
    ${SWEEP_ARGS[@]+"${SWEEP_ARGS[@]}"} "$@")
merge_jid=$(bicamera_submit_array "${LABEL}_merge" 1 "${sweep_jid}" "" \
    benchmarks/merge_shards.sh "${SHARD_DIR}" "${LABEL}" "${N}" "${MERGED}")
popd >/dev/null

cat <<EOF
BICAMERA benchmark sweep
  sweep:        ${SWEEP}$( [ "${SWEEP_ARGS[0]:-}" = "--grid" ] && echo " (grid ${GRID})" )
  script:       benchmarks/${SCRIPT}
  tasks:        ${N}
  extra args:   ${*:-none}
  shards:       ${SHARD_DIR}/${LABEL}_task_<i>.tsv
  merged:       ${MERGED}
  scheduler:    ${BICAMERA_SCHEDULER}
EOF
bicamera_scheduler_epilogue
