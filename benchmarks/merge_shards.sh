#!/bin/bash
# ABOUTME: Merges a sweep's per-task shards into one TSV: header once, then shards 1..N in task order.
# ABOUTME: Queued by run_benchmarks.sh after the array. Args: <shard_dir> <label> <n_tasks> <out_tsv>.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
set -euo pipefail

SHARD_DIR=$1
LABEL=$2
N=$3
OUT=$4

# Refuse to merge a partial sweep: a missing shard means a failed task, and a
# table without it would silently drop a whole depth (or n_asm) from every figure.
missing=0
for i in $(seq 1 "${N}"); do
    f="${SHARD_DIR}/${LABEL}_task_${i}.tsv"
    [ -s "${f}" ] || { echo "ERROR: missing or empty shard ${f}" >&2; missing=1; }
done
[ "${missing}" -eq 0 ] || exit 1

tmp="${OUT}.tmp"
head -1 "${SHARD_DIR}/${LABEL}_task_1.tsv" > "${tmp}"
for i in $(seq 1 "${N}"); do
    tail -n +2 "${SHARD_DIR}/${LABEL}_task_${i}.tsv" >> "${tmp}"
done
mv "${tmp}" "${OUT}"

echo "Merged ${N} ${LABEL} shards -> ${OUT} ($(($(wc -l < "${OUT}") - 1)) rows)"
