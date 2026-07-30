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

export CHROMS="${CHROMS:-$(seq -f "chr%g" 1 22) chrX chrY}"
export WINDOWSIZE="${WINDOWSIZE:-10}"
export WINDOW_STEP="${WINDOW_STEP:-5}"
export WORKERS="${WORKERS:-10}"
# ASM restart override. Empty default → driver picks its built-in default (5).
export ASM_RESTARTS="${ASM_RESTARTS:-}"
# ASM (M1) EM convergence tol override. Empty default → driver/library default (1e-4).
export EM_TOL="${EM_TOL:-}"
# Potential-ASM-region front-end tunables (pileup -> select -> merge).
# A CpG is a candidate iff valid coverage >= MIN_COV and FRAC_LOW <= population
# methylation fraction <= FRAC_HIGH; candidates within MERGE_DIST bp merge into
# one potential-ASM region. These feed asm_regions_inner.sh via qsub -V.
export MIN_COV="${MIN_COV:-10}"
# FRAC_LOW/FRAC_HIGH bound the pooled (phase-ignored) methylation fraction. The
# defaults are the population-mean envelope of a balanced hypo/hyper site under the
# Atlas (Nat Commun 2025) per-fragment cutoffs (hypo <=0.35, hyper >=0.65):
#   0.5*0.65 = 0.325 (lower)   0.5*(1+0.35) = 0.675 (upper).
# See docs/cpg_selection_threshold_guide.md for the full rationale.
export FRAC_LOW="${FRAC_LOW:-0.325}"
export FRAC_HIGH="${FRAC_HIGH:-0.675}"
export MERGE_DIST="${MERGE_DIST:-1000}"
# Drop potential-ASM regions with fewer than this many distinct in-region CpGs
# (bucket step). Regions with >= this but < WINDOWSIZE get one whole-region EM window.
export MIN_REGION_CPGS="${MIN_REGION_CPGS:-3}"
RUN_LABEL="${RUN_LABEL:-bam}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Repo root, exported so the inner jobs can locate scripts/ and em/. The inners
# cannot derive it themselves: under SGE, $0 / BASH_SOURCE[0] resolve to
# /var/spool/ge/.../job_scripts/<JID>. They previously used SGE_O_WORKDIR, which
# only exists under SGE and so breaks the local scheduler backend.
export BICAMERA_ROOT="${SCRIPT_DIR}"
export ASM_DRIVER="${BICAMERA_ROOT}/em/BIC_algorithm_ASM.py"
export ASM_OUT="${OUT_DIR}/ASM"

# modkit binary — exported so SGE inner scripts resolve correctly
# (under SGE, $0 / BASH_SOURCE[0] resolve to /var/spool/ge/.../job_scripts/).
export MODKIT_BIN="${SCRIPT_DIR}/bin/modkit"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------
for f in "${PHASED_BAM}" "${REF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: file not found: ${f}" >&2
        exit 1
    fi
done
if [ ! -x "${MODKIT_BIN}" ]; then
    echo "ERROR: modkit binary not executable: ${MODKIT_BIN}" >&2
    exit 1
fi
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
    if command -v samtools >/dev/null 2>&1; then
        # Locally disable pipefail: head closes the pipe after 5000 reads, which
        # SIGPIPE-kills samtools. We only need to know if MM was seen.
        set +o pipefail
        found=$(samtools view -F 260 "${PHASED_BAM}" 2>/dev/null | head -5000 | awk '
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
    else
        echo "WARN: samtools not on PATH — skipping modBAM tag check. Load samtools/1.19 before submission to enable." >&2
    fi
fi

# ------------------------------------------------------------
# Set up sample working dir with phased BAM symlinks
# ------------------------------------------------------------
mkdir -p "${OUT_DIR}" "${OUT_DIR}/log" "${ASM_OUT}" "${SCRIPT_DIR}/log"
ln -sfn "$(readlink -f "${PHASED_BAM}")" "${OUT_DIR}/phased_output.bam"
if [ -f "${PHASED_BAM}.bai" ]; then
    ln -sfn "$(readlink -f "${PHASED_BAM}.bai")" "${OUT_DIR}/phased_output.bam.bai"
elif [ -f "${PHASED_BAM%.bam}.bai" ]; then
    ln -sfn "$(readlink -f "${PHASED_BAM%.bam}.bai")" "${OUT_DIR}/phased_output.bam.bai"
else
    echo "WARN: no .bai index alongside ${PHASED_BAM}" >&2
fi

# ------------------------------------------------------------
# Submit the whole DAG as SGE array jobs up front — one task per chromosome per
# stage — chained by captured job ID. Two legs run in parallel and join at modkitref:
#   discovery : pileup -> regions      (potential-ASM region BEDs)
#   read      : hpmap                  (read_id->HP map straight from the original BAM)
# Both hpmap and modkitref read the indexed original phased_output.bam with --region
# <chr> (the per-chr split BAM stage was deleted — see CHANGES.md 2026-06-30), so
# hpmap is now a leaf that can start immediately on the launch-time BAM symlink.
# modkitref holds on BOTH hpmap and regions (read-level extraction is restricted to
# the candidate regions); then: modkitref -> filter -> coordfix -> m2r -> bucket -> ASM EM.
#
# Every inner derives its chromosome from SGE_TASK_ID indexing ${CHROMS}, so every
# submit MUST pass -V (CHROMS is exported above) and the array is sized -t 1-N from
# the SAME list. Job IDs are captured with -terse (array form prints
# <jobid>.<start>-<end>:<step>; cut keeps the bare id); each stage holds on the
# previous stage's id. ID holds are run-scoped (immune to RUN_LABEL name collisions)
# and, unlike a name hold, do not release on a crashed predecessor task.
# Per-inner resources (s_vmem, def_slot) stay in each inner's #$ directives.
# ------------------------------------------------------------
chromosomes=(${CHROMS})
N=${#chromosomes[@]}
# Per-slot memory for the ASM EM array (default 2G; raise via ASM_EM_VMEM=4G for
# ONT/unfiltered inputs — see memory_rightsizing_report.md).
ASM_EM_VMEM="${ASM_EM_VMEM:-2G}"
bare_id() { cut -d. -f1; }   # bare job id from a -terse array submit

pushd "${SCRIPT_DIR}" >/dev/null

# discovery leg: pileup -> regions
pileup_jid=$(qsub -terse -V -N "${RUN_LABEL}_pileup" -t 1-${N} \
    scripts/asm_pileup_inner.sh "${OUT_DIR}" "${REF}" | bare_id)
regions_jid=$(qsub -terse -V -N "${RUN_LABEL}_regions" -t 1-${N} -hold_jid "${pileup_jid}" \
    scripts/asm_regions_inner.sh "${OUT_DIR}" | bare_id)

# read leg: hpmap (reads the indexed original BAM directly with a region arg; no split).
# Leaf job — only needs the launch-time phased_output.bam(.bai) symlink, so no hold.
hpmap_jid=$(qsub -terse -V -N "${RUN_LABEL}_hpmap" -t 1-${N} \
    scripts/hpmap_inner.sh "${OUT_DIR}" | bare_id)

# join + shared chain: modkitref -> filter -> coordfix -> m2r -> bucket -> ASM EM
modkitref_jid=$(qsub -terse -V -N "${RUN_LABEL}_modkitref" -t 1-${N} \
    -hold_jid "${hpmap_jid},${regions_jid}" \
    scripts/modkitinner_withreference_bychrandhap.sh "${OUT_DIR}" "${REF}" | bare_id)
filter_jid=$(qsub -terse -V -N "${RUN_LABEL}_filter" -t 1-${N} -hold_jid "${modkitref_jid}" \
    scripts/modkitfilterinner_bychrandhap.sh "${OUT_DIR}" | bare_id)
fixer_jid=$(qsub -terse -V -N "${RUN_LABEL}_fixer" -t 1-${N} -hold_jid "${filter_jid}" \
    scripts/coordfixer_inner.sh "${OUT_DIR}" | bare_id)
m2r_jid=$(qsub -terse -V -N "${RUN_LABEL}_m2r" -t 1-${N} -hold_jid "${fixer_jid}" \
    scripts/merged2reads_inner.sh "${OUT_DIR}" | bare_id)
bucket_jid=$(qsub -terse -V -N "${RUN_LABEL}_bucket" -t 1-${N} -hold_jid "${m2r_jid}" \
    scripts/asmregion_bucket_inner.sh "${OUT_DIR}" | bare_id)
asm_jid=$(qsub -terse -V -N "${RUN_LABEL}_asm" -t 1-${N} -hold_jid "${bucket_jid}" \
    -l s_vmem=${ASM_EM_VMEM},mem_req=${ASM_EM_VMEM} \
    scripts/run_ASM_inner.sh "${ASM_OUT}" "${OUT_DIR}" | bare_id)

popd >/dev/null

cat <<EOF
Submitted ASM pipeline (from BAM, potential-ASM-region route).
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

Monitor:  qstat -u \$USER | grep ${RUN_LABEL}
EOF

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
