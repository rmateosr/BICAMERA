#!/bin/bash
# ABOUTME: Build a per-chr read_id->HP map (H1/H2/noH) from PRIMARY alignments only.
# ABOUTME: Replaces the 3-way BAM split; the map is joined onto the single modkit extract at coordfix.
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -l s_vmem=2G
set -xv
set -o errexit
set -o nounset

module use /usr/local/package/modulefiles/
module load samtools/1.19

OUTPUT_FOLDER=$1
# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
chr=${chromosomes[$((SGE_TASK_ID-1))]}
mkdir -p ${OUTPUT_FOLDER}/hpmap

# modkit extract calls emits rows for PRIMARY alignments only, so a primary-only
# map (-F 0x900 excludes secondary 0x100 + supplementary 0x800) covers exactly the
# read_ids that appear in the extract and is guaranteed unique (one primary per read).
# HP:i:1 -> H1, HP:i:2 -> H2, anything else (no HP tag) -> noH. This reproduces the
# old file-provenance partition exactly (verified row-set identical on chr15 HG002).
#
# Read the indexed original BAM with a region arg (chr) instead of a per-chr split
# BAM: `samtools view -F 0x900 <bam> <chr>` index-seeks to chr and yields exactly the
# primary records the (now-deleted) splitbam stage produced, so the read_id->HP map is
# unchanged. The .bai symlink alongside phased_output.bam is made by run_pipeline.sh.
OUT=${OUTPUT_FOLDER}/hpmap/${chr}.tsv
samtools view -F 0x900 ${OUTPUT_FOLDER}/phased_output.bam ${chr} | awk -F'\t' '
{
    hp="noH";
    for(i=12;i<=NF;i++){ if($i ~ /^HP:i:/){ v=substr($i,6); if(v=="1")hp="H1"; else if(v=="2")hp="H2"; } }
    print $1"\t"hp;
}' > ${OUT}.tmp
mv ${OUT}.tmp ${OUT}
