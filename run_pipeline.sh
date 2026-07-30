#!/bin/bash
# ABOUTME: Static ASM pipeline starting from a phased modBAM (modkit path), single sample.
# ABOUTME: pileup->select->merge regions; modkit(--region+--include-bed)->filter->coordfix->m2r->bucket->Shared+ASM BIC.
set -euo pipefail

# ============================================================
# Usage
# ============================================================
# bash run_pipeline.sh <PHASED_BAM> <REF> <OUT_DIR>
#
# PHASED_BAM : Clair3+WhatsHap output with HP tags and MM/ML modification tags.
# REF        : reference genome FASTA (required by modkit and downstream scripts).
# OUT_DIR    : per-sample working directory (created if missing).
#
# Optional env overrides:
#   CHROMS              — space-separated chromosomes (default: chr1..22 chrX chrY)
#   WINDOWSIZE          — CpG sites per EM window (default: 10)
#   WINDOW_STEP         — CpGs advanced between EM windows (default: 5; 1 = every CpG)
#   WORKERS             — parallel processes for ASM EM driver (default: 10)
#   RUN_LABEL           — prefix for SGE job names (default: bam)
#   ASM_RESTARTS        — EM restarts for the ASM BIC; unset = driver default 5
#   EM_TOL              — ASM (M1) EM convergence tol; unset = driver default 1e-4
#
# NOTE: this pipeline requires modBAM input (Dorado-style MM/ML tags).
# Guppy/nanopolish-era basecalls are not supported. The startup guard below
# aborts if MM/ML are absent.
# ============================================================

if [ "$#" -ne 3 ]; then
    cat <<EOF
Usage: $0 <PHASED_BAM> <REF> <OUT_DIR>
EOF
    exit 2
fi

PHASED_BAM=$1
REF=$2
OUT_DIR=$3

# Resolve to absolute paths. The inner jobs run with cwd = the repo root, so a
# RELATIVE OUT_DIR/REF passed here would resolve against the repo inside every
# job and silently break with "file not found". realpath -m does not require the
# path to exist yet.
REF="$(realpath -m "${REF}")"
OUT_DIR="$(realpath -m "${OUT_DIR}")"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Repo root, exported so the inner jobs can locate scripts/, config/ and em/. The
# inners cannot derive it themselves: under SGE, $0 / BASH_SOURCE[0] resolve to
# /var/spool/ge/.../job_scripts/<JID>. They previously used SGE_O_WORKDIR, which
# only exists under SGE and so breaks the local scheduler backend.
export BICAMERA_ROOT="${SCRIPT_DIR}"

# Chromosomes to process. Defined HERE and nowhere else: this list sizes the job
# array (-t 1-N) and each inner derives its own chromosome by indexing the SAME
# list, so a second independently-defaulted definition could make an inner
# process a different chromosome than the one its task index was meant for.
# The default is UCSC chr-prefixed human; override for other assemblies or
# naming conventions, e.g. CHROMS="1 2 3" for Ensembl-style names.
export CHROMS="${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}"

# Tool resolution (modkit/samtools/bedtools/python) and every tunable default,
# both shared verbatim with the inner jobs so there is exactly one definition of
# each. See config/site.sh and config/defaults.sh.
. "${BICAMERA_ROOT}/config/site.sh"
. "${BICAMERA_ROOT}/config/defaults.sh"

# Scheduler backend. Defaults to sge where qsub exists and local otherwise, so the
# pipeline runs on a workstation with no configuration. Force one with
# BICAMERA_SCHEDULER=sge|local; add a backend by dropping a scheduler/<name>.sh
# that defines bicamera_submit_array and bicamera_scheduler_epilogue.
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
    echo "       available: $(cd "${BICAMERA_ROOT}/scheduler" && ls *.sh | sed 's/\.sh$//' | tr '\n' ' ')" >&2
    exit 1
fi
. "${_bicamera_backend}"
unset _bicamera_backend

export ASM_DRIVER="${BICAMERA_ROOT}/em/BIC_algorithm_ASM.py"
export ASM_OUT="${OUT_DIR}/ASM"
RUN_LABEL="${RUN_LABEL:-bam}"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------
for f in "${PHASED_BAM}" "${REF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: file not found: ${f}" >&2
        exit 1
    fi
done
# Check every tool the DAG will need up front, so a missing bedtools surfaces now
# rather than after the 24 pileup tasks have already run.
bicamera_require modkit samtools bedtools python
if [ ! -f "${ASM_DRIVER}" ]; then
    echo "ERROR: ASM driver not found: ${ASM_DRIVER}" >&2
    exit 1
fi

# ------------------------------------------------------------
# modBAM guard: ensure inputs carry MM/ML methylation tags.
# Empty-tag BAMs (Guppy/nanopolish-era) silently produce empty modkit TSVs
# and crash coordfixer mid-chain. Sample the first 5000 mapped reads;
# abort with a clear error if no MM tag is present.
# Set SKIP_MODBAM_CHECK=1 to bypass entirely.
# ------------------------------------------------------------
if [ "${SKIP_MODBAM_CHECK:-0}" != "1" ]; then
        # Locally disable pipefail: head closes the pipe after 5000 reads, which
        # SIGPIPE-kills samtools. We only need to know if MM was seen.
        set +o pipefail
        found=$("${SAMTOOLS}" view -F 260 "${PHASED_BAM}" 2>/dev/null | head -5000 | awk '
            { for(i=12;i<=NF;i++) if($i ~ /^M[Mm]:/) { print 1; exit } }
        ')
        set -o pipefail
        if [ "${found:-0}" != "1" ]; then
            cat <<ERR >&2
ERROR: Input BAM has no MM/ML methylation tags in the first 5000 reads.
       file : ${PHASED_BAM}
       BICAMERA requires modBAM input carrying Dorado-style MM/ML tags.
       Guppy / nanopolish-era basecalls are not supported.
       To bypass this check, set SKIP_MODBAM_CHECK=1.
ERR
            exit 3
        fi
fi

# ------------------------------------------------------------
# Set up sample working dir with phased BAM symlinks
# ------------------------------------------------------------
# Logs go under OUT_DIR, not into the repo — the backend passes -o/-e there.
mkdir -p "${OUT_DIR}" "${OUT_DIR}/log" "${ASM_OUT}"
ln -sfn "$(readlink -f "${PHASED_BAM}")" "${OUT_DIR}/phased_output.bam"
if [ -f "${PHASED_BAM}.bai" ]; then
    ln -sfn "$(readlink -f "${PHASED_BAM}.bai")" "${OUT_DIR}/phased_output.bam.bai"
elif [ -f "${PHASED_BAM%.bam}.bai" ]; then
    ln -sfn "$(readlink -f "${PHASED_BAM%.bam}.bai")" "${OUT_DIR}/phased_output.bam.bai"
else
    echo "WARN: no .bai index alongside ${PHASED_BAM}" >&2
fi

# ------------------------------------------------------------
# Run the whole DAG — one task per chromosome per stage. Two legs are independent
# and join at modkitref:
#   discovery : pileup -> regions      (potential-ASM region BEDs)
#   read      : hpmap                  (read_id->HP map straight from the original BAM)
# Both hpmap and modkitref read the indexed original phased_output.bam with
# --region <chr> (there is no per-chr split BAM stage), so hpmap is a leaf that can
# start immediately on the launch-time BAM symlink. modkitref waits on BOTH hpmap
# and regions (read-level extraction is restricted to the candidate regions); then
# modkitref -> filter -> coordfix -> m2r -> bucket -> ASM EM.
#
# Stage submission goes through bicamera_submit_array, provided by the selected
# scheduler backend. The SGE backend queues every stage up front and chains them by
# job id, so this launcher returns immediately; the local backend runs each stage to
# completion in call order, which satisfies the same dependencies serially. The
# dependency and resource arguments are therefore advisory to the backend — see
# scheduler/local.sh for what it ignores and why.
#
# Every inner derives its chromosome from SGE_TASK_ID indexing ${CHROMS}, so each
# array is sized from the SAME list the inners will index.
# ------------------------------------------------------------
chromosomes=(${CHROMS})
N=${#chromosomes[@]}

pushd "${BICAMERA_ROOT}" >/dev/null

# discovery leg: pileup -> regions
pileup_jid=$(bicamera_submit_array pileup "${N}" "" "" \
    scripts/asm_pileup_inner.sh "${OUT_DIR}" "${REF}")
regions_jid=$(bicamera_submit_array regions "${N}" "${pileup_jid}" "" \
    scripts/asm_regions_inner.sh "${OUT_DIR}")

# read leg: hpmap (reads the indexed original BAM directly with a region arg; no split).
# Leaf stage — only needs the launch-time phased_output.bam(.bai) symlink.
hpmap_jid=$(bicamera_submit_array hpmap "${N}" "" "" \
    scripts/hpmap_inner.sh "${OUT_DIR}")

# join + shared chain: modkitref -> filter -> coordfix -> m2r -> bucket -> ASM EM
modkitref_jid=$(bicamera_submit_array modkitref "${N}" "${hpmap_jid},${regions_jid}" "" \
    scripts/modkitinner_withreference_bychrandhap.sh "${OUT_DIR}" "${REF}")
filter_jid=$(bicamera_submit_array filter "${N}" "${modkitref_jid}" "" \
    scripts/modkitfilterinner_bychrandhap.sh "${OUT_DIR}")
fixer_jid=$(bicamera_submit_array fixer "${N}" "${filter_jid}" "" \
    scripts/coordfixer_inner.sh "${OUT_DIR}")
m2r_jid=$(bicamera_submit_array m2r "${N}" "${fixer_jid}" "" \
    scripts/merged2reads_inner.sh "${OUT_DIR}")
bucket_jid=$(bicamera_submit_array bucket "${N}" "${m2r_jid}" "" \
    scripts/asmregion_bucket_inner.sh "${OUT_DIR}")
asm_jid=$(bicamera_submit_array asm "${N}" "${bucket_jid}" \
    "s_vmem=${ASM_EM_VMEM},mem_req=${ASM_EM_VMEM}" \
    scripts/run_ASM_inner.sh "${ASM_OUT}" "${OUT_DIR}")

popd >/dev/null

cat <<EOF
BICAMERA — static ASM detection from a phased modBAM.
  BAM:            ${PHASED_BAM}
  reference:      ${REF}
  modkit:         ${MODKIT_BIN}
  output dir:     ${OUT_DIR}
  regions BED:    ${OUT_DIR}/potential_asm_regions/<chr>.bed
  ASM TSVs:       ${ASM_OUT}/<chr>_ASM.tsv
  ASM driver:     ${ASM_DRIVER}
  chromosomes:    ${CHROMS}
  windowsize:     ${WINDOWSIZE}
  window step:    ${WINDOW_STEP}
  workers:        ${WORKERS}
  asm restarts:   ${ASM_RESTARTS:-driver default}
  em tol:         ${EM_TOL:-driver default}
  region select:  cov>=${MIN_COV}, ${FRAC_LOW}<=frac<=${FRAC_HIGH}, merge -d ${MERGE_DIST}
  run label:      ${RUN_LABEL}
  site profile:   ${BICAMERA_SITE}
  scheduler:      ${BICAMERA_SCHEDULER}
EOF
bicamera_scheduler_epilogue

# ------------------------------------------------------------
# Append a run-registry entry so we can track which runs have been done,
# when, and at what code level. Records date, sample BAM, reference, output
# dir, tunables, the reproducible command, and the git commit (+ whether the
# pipeline code differs from HEAD). Best-effort: a logging hiccup must never
# fail an already-submitted run, so the whole block falls through to a WARN.
# ------------------------------------------------------------
RUN_LOG="${SCRIPT_DIR}/RUN_LOG.md"
{
    RUN_TS="$(date '+%Y-%m-%d %H:%M:%S %Z')"
    # CHROMS default is built with `seq` (one per line); collapse to a single line.
    CHROMS_FMT="$(echo ${CHROMS})"
    if GIT_SHA="$(git -C "${SCRIPT_DIR}" rev-parse --short HEAD 2>/dev/null)"; then
        GIT_BRANCH="$(git -C "${SCRIPT_DIR}" rev-parse --abbrev-ref HEAD 2>/dev/null)"
        # "dirty" only w.r.t. TRACKED pipeline code. Use -uno so untracked files
        # and the run log itself don't trip the flag — it means "committed pipeline
        # code differs from HEAD", which is what "are these runs at the same level"
        # needs. SCRIPT_DIR is the repo root, so this covers em/ as well.
        if [ -n "$(git -C "${SCRIPT_DIR}" status --porcelain -uno 2>/dev/null | grep -v 'RUN_LOG.md')" ]; then
            GIT_STATE="${GIT_SHA} (${GIT_BRANCH}, pipeline code DIRTY vs HEAD)"
        else
            GIT_STATE="${GIT_SHA} (${GIT_BRANCH}, clean)"
        fi
    else
        GIT_STATE="unknown (not a git checkout)"
    fi
    if [ ! -f "${RUN_LOG}" ]; then
        cat > "${RUN_LOG}" <<'HDR'
# RUN_LOG — BICAMERA pipeline runs

Auto-appended by `run_pipeline.sh` on every submission (newest at the bottom).
Each entry records the sample BAM, reference, output dir, tunables, the exact
reproducible command, and the git commit — so runs can be confirmed to share a
code "level". A `DIRTY vs HEAD` git state means the pipeline code had uncommitted
changes at run time, so that run is NOT pinned to the commit.

---
HDR
    fi
    cat >> "${RUN_LOG}" <<EOF

## ${RUN_TS} — ${RUN_LABEL}

- **Sample BAM:** \`${PHASED_BAM}\`
- **Reference:** \`${REF}\`
- **Output dir:** \`${OUT_DIR}\`
- **Chromosomes:** ${CHROMS_FMT}
- **Git commit:** ${GIT_STATE}
- **Tunables:** windowsize=${WINDOWSIZE} step=${WINDOW_STEP} workers=${WORKERS} restarts=${ASM_RESTARTS:-default} tol=${EM_TOL:-default}; region cov>=${MIN_COV} ${FRAC_LOW}<=frac<=${FRAC_HIGH} merge=${MERGE_DIST} min_region_cpgs=${MIN_REGION_CPGS}
- **Command:** \`RUN_LABEL=${RUN_LABEL} bash ${SCRIPT_DIR}/run_pipeline.sh "${PHASED_BAM}" "${REF}" "${OUT_DIR}"\`
EOF
    echo "Logged run to: ${RUN_LOG}"
} || echo "WARN: failed to append run-registry entry to ${RUN_LOG}" >&2
