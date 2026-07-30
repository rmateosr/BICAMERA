#!/bin/bash
# ABOUTME: Per-chr read-level modkit extract, restricted to the potential-ASM regions.
# ABOUTME: One pass over the indexed BAM (--region <chr> + --include-bed) -> per-read CpG calls.
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=12G
set -xv
set -o errexit
set -o nounset

# Repo root: exported by run_pipeline.sh. SGE_O_WORKDIR / $(pwd) cover a stage
# invoked by hand. site.sh resolves the tools; defaults.sh supplies the tunables.
BICAMERA_ROOT="${BICAMERA_ROOT:-${SGE_O_WORKDIR:-$(pwd)}}"
. "${BICAMERA_ROOT}/config/site.sh"
. "${BICAMERA_ROOT}/config/defaults.sh"
bicamera_require modkit

OUTPUT_FOLDER=$1
REF=$2
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
chr=${chromosomes[$((SGE_TASK_ID-1))]}
# Output dir was created by the (now-removed) outer; create it here so the task is self-contained.
mkdir -p "${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap"
# modkit v0.4.1 moved read-level extraction under the `calls` subcommand.
# Top-level `modkit extract` no longer accepts flags directly — it requires a
# subcommand (full | calls). `calls` preserves the --cpg / --mapped-only /
# --reference flags and emits a read-level TSV with the same column layout
# downstream filters expect ($3 = ref_position, $12 = call_code).
# --threads 1 here is a REPRODUCIBILITY choice, not slot-matching. Threading
# `modkit extract calls` is a PURE REORDER of the output: each read's CpG rows still
# come out as one contiguous block (verified 0 non-contiguous reads at 4 AND 8 threads,
# 8k reads), so per-read aggregation is NOT corrupted. But the BETWEEN-read block order
# follows thread completion instead of BAM order, and that read-order propagates to the
# EM matrix X row order. The M1 EM seed is md5(X.tobytes()) (see EMfunctions_ASM.py
# _seed_from_X), so a different read order -> different seed -> different EM trajectory
# -> different BIC vs the committed --threads 1 baseline (verified 2026-06-30: threaded
# chr15 ASM has a different sorted-md5; M0/null is order-invariant, only M1 shifts).
# Threading is therefore safe ONLY paired with a canonical read-order sort + a one-time
# re-baseline, or by making the EM seed order-invariant (an asm-stage change). Until
# then keep --threads 1. (pileup IS threaded — its per-site counts are order-independent.)
#
# --no-filtering: keep EVERY call and binarize at a fixed 0.5 (the call_code is
# already modkit's 0.5 argmax: m iff P(5mC)>0.5, - otherwise). We deliberately do
# NOT use modkit's default adaptive threshold (the 10th-percentile pass/fail), which
# is estimated by sampling the BAM and therefore depends on which reads are present
# (per-haplotype / per-condition drift). A fixed 0.5 is composition-independent and
# reproducible, and matches the ASMS reference tool (binarize at 127/255 ≈ 0.5, no
# confidence drop). The downstream filter no longer reads the `fail` column.
#
# Potential-ASM-region route: restrict read-level extraction to the
# methylation-guided regions (potential_asm_regions/<chr>.bed).
#   * BED present and non-empty -> --include-bed restricts extraction.
#   * BED present but EMPTY (chromosome with no intermediate-methylation CpGs) ->
#     valid state; write an empty extract and skip modkit (it errors on an empty
#     --include-bed).
#   * BED absent (route run without the region front-end) -> whole-chromosome.
# Single pass over the indexed original BAM, restricted to this chromosome with
# --region ${chr} (no per-chr split BAM). --region (coarse read-fetch by alignment
# position) composes with --include-bed (fine per-CpG site filter); verified
# byte-identical to the old split-BAM + --include-bed extract on a chr15 region.
# modkit extract calls only emits primary alignments, so the per-row read_id matches
# the primary-only hpmap built upstream; the haplotype label is joined in at coordfix
# instead of encoded by file provenance.
REGION_BED="${OUTPUT_FOLDER}/potential_asm_regions/${chr}.bed"
OUT_TSV="${OUTPUT_FOLDER}/modkit_splitbam_referenced_bychrandhap/modkit_extract_${chr}.tsv"
if [ -f "${REGION_BED}" ] && [ ! -s "${REGION_BED}" ]; then
    echo "NOTE: empty region BED for ${chr}; writing empty extract ${OUT_TSV}." >&2
    : > "${OUT_TSV}"
else
    INCLUDE_FLAGS=()
    if [ -s "${REGION_BED}" ]; then
        INCLUDE_FLAGS+=(--include-bed "${REGION_BED}")
    fi
    "${MODKIT_BIN}" extract calls --threads 1 --cpg --mapped-only --no-filtering --region "${chr}" "${INCLUDE_FLAGS[@]}" --reference "${REF}" "${OUTPUT_FOLDER}/phased_output.bam" "${OUT_TSV}"
fi




