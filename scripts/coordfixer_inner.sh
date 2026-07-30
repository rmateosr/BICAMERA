#!/bin/bash
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

# Array task -> chromosome: index the SAME CHROMS list used to size -t 1-N (via -V).
: "${CHROMS:?CHROMS not in task env — submit the array with -V}"
chromosomes=(${CHROMS})
CHR=${chromosomes[$((SGE_TASK_ID-1))]}
echo ${CHR}	
OUTPUT_FOLDER=$1
echo ${OUTPUT_FOLDER}	
python3.12 scripts/coordfixer.py ${CHR} ${OUTPUT_FOLDER}



