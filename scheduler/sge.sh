# ABOUTME: SGE backend — submits each pipeline stage as a per-chromosome array job.
# ABOUTME: Stages are chained by captured job id, so the whole DAG is queued up front.

# bicamera_submit_array <stage> <n_tasks> <hold_handles> <resources> <script> [args...]
#
# Queues <script> as an array of <n_tasks> tasks and prints the job id. Each task
# derives its chromosome from SGE_TASK_ID by indexing $CHROMS, so the array MUST be
# sized from the same list the inners will index, and -V is mandatory (it carries
# CHROMS, BICAMERA_ROOT and the tunables).
#
#   <hold_handles>  comma-separated job ids to wait on, or "" for none. ID holds are
#                   run-scoped, so they are immune to RUN_LABEL collisions between
#                   concurrent runs and — unlike name holds — do not release when a
#                   predecessor task crashes.
#   <resources>     extra `-l` spec (e.g. "s_vmem=4G,mem_req=4G"), or "" to use only
#                   the #$ directives in the inner. These complex names are
#                   site-defined; the local backend ignores them.
#
# -terse on an array submit prints "<jobid>.<start>-<end>:<step>", so cut keeps the
# bare id. Logs are directed to $OUT_DIR/log rather than a path inside the repo.
bicamera_submit_array() {
    local stage="$1" n="$2" holds="$3" resources="$4"
    shift 4
    local qsub_args=(
        -terse -V
        -N "${RUN_LABEL}_${stage}"
        -t "1-${n}"
        -o "${OUT_DIR}/log" -e "${OUT_DIR}/log"
    )
    [ -n "${holds}" ]     && qsub_args+=(-hold_jid "${holds}")
    [ -n "${resources}" ] && qsub_args+=(-l "${resources}")
    qsub "${qsub_args[@]}" "$@" | cut -d. -f1
}

# Printed once after all stages are queued.
bicamera_scheduler_epilogue() {
    cat <<EOF

Backend: SGE (jobs queued; this launcher returns immediately).
Monitor:  qstat -u \$USER | grep ${RUN_LABEL}
Logs:     ${OUT_DIR}/log
EOF
}
