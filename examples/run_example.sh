#!/bin/bash
# ABOUTME: Worked example — run BICAMERA on your own phased modBAM, one chromosome or all.
# ABOUTME: Copy and edit, or just pass your paths as arguments.
set -euo pipefail

# ============================================================
# Usage
# ============================================================
#   bash examples/run_example.sh <PHASED_BAM> <REF> <OUT_DIR> [CHR ...]
#
# PHASED_BAM  modBAM with Dorado-style MM/ML methylation tags and HP haplotype
#             tags (e.g. Clair3 + WhatsHap output). Must be indexed (.bai
#             alongside). See docs/input_requirements.md for how to produce one.
# REF         reference FASTA the BAM was aligned to.
# OUT_DIR     working directory for this run; created if missing. Everything the
#             run produces lands here, including its RUN_LOG.md. Keep it OUTSIDE
#             the repo — intermediates reach ~130 GB per chromosome.
# CHR ...     optional chromosome list. Default: chr1..chr22, chrX, chrY.
#
# Start with ONE small chromosome to confirm the setup before a whole genome:
#   bash examples/run_example.sh sample.bam ref.fa /scratch/asm_chr22 chr22
#
# Tunables are environment variables; every default lives in config/defaults.sh.
# The ones worth knowing:
#   WINDOWSIZE=10       CpG sites per EM window
#   WINDOW_STEP=5       CpGs advanced between windows (1 = every CpG)
#   WORKERS=10          parallel EM processes
#   MIN_COV=10          minimum coverage for a candidate CpG
#   FRAC_LOW/FRAC_HIGH  pooled-methylation band for candidate CpGs (0.325/0.675)
#   MERGE_DIST=1000     bp within which candidate CpGs merge into one region
#
# Environment:
#   BICAMERA_SCHEDULER=sge|local   default: sge if qsub exists, else local
#   BICAMERA_SITE=<name>|none      load config/site.<name>.sh; none = PATH only
#   RUN_LABEL=<prefix>             job-name prefix and log prefix
# ============================================================

if [ "$#" -lt 3 ]; then
    sed -n '6,38p' "${BASH_SOURCE[0]}" | sed 's/^# \?//' >&2
    exit 2
fi

PHASED_BAM=$1
REF=$2
OUT_DIR=$3
shift 3

REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

# Any remaining arguments are the chromosome list. Leaving CHROMS unset lets
# run_pipeline.sh apply its default; note it must be a single space-separated
# string, because that string both sizes the job array and is indexed per task.
if [ "$#" -gt 0 ]; then
    export CHROMS="$*"
fi

export RUN_LABEL="${RUN_LABEL:-example}"

exec bash "${REPO_ROOT}/run_pipeline.sh" "${PHASED_BAM}" "${REF}" "${OUT_DIR}"
