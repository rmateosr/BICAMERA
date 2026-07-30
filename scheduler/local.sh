# ABOUTME: Local backend — runs each pipeline stage serially, no scheduler needed.
# ABOUTME: Same interface as scheduler/sge.sh; for workstations and small test runs.

# bicamera_submit_array <stage> <n_tasks> <hold_handles> <resources> <script> [args...]
#
# Runs <script> once per task, in order, setting SGE_TASK_ID so the inners derive
# their chromosome exactly as they do under SGE. Blocks until the stage finishes —
# that is what makes the dependency arguments unnecessary here:
#
#   <hold_handles>  ignored. run_pipeline.sh calls the stages in dependency order,
#                   and every call blocks, so the ordering the holds encode is
#                   already satisfied. Note this serializes the two parallel legs
#                   (pileup->regions and hpmap) — slower than SGE, same results.
#   <resources>     ignored. s_vmem / mem_req are SGE complexes with no local
#                   equivalent; a local run is limited by the machine.
#
# Per-task output goes to $OUT_DIR/log/<label>_<stage>.o<task>, mirroring the SGE
# layout. A failing task aborts the stage immediately and names its log, rather
# than letting later stages run on incomplete input.
bicamera_submit_array() {
    local stage="$1" n="$2" holds="$3" resources="$4"
    shift 4
    local script="$1"
    shift
    local i log rc
    mkdir -p "${OUT_DIR}/log"
    for i in $(seq 1 "${n}"); do
        log="${OUT_DIR}/log/${RUN_LABEL}_${stage}.o${i}"
        # Capture the task's status with `|| rc=$?`, NOT `if ! cmd`: inside an
        # `if !` the value of $? is the negation's status, i.e. always 0, so a
        # failing task would report exit 0 and this function would return success
        # while the launcher carried on to the next stage.
        rc=0
        SGE_TASK_ID="${i}" bash "${script}" "$@" >"${log}" 2>&1 || rc=$?
        if [ "${rc}" -ne 0 ]; then
            echo "ERROR: stage '${stage}' task ${i}/${n} failed (exit ${rc})." >&2
            echo "       log: ${log}" >&2
            return "${rc}"
        fi
    done
    # The launcher captures this as a handle; nothing downstream needs to wait.
    echo "local:${stage}"
}

bicamera_scheduler_epilogue() {
    cat <<EOF

Backend: local (all stages already ran to completion).
Logs:     ${OUT_DIR}/log
EOF
}
