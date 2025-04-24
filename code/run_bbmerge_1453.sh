#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --account=gbru_fy24_stinkbug_diet
#SBATCH --job-name=bbmerge_1453
#SBATCH --output=bbmerge_1453_%j.out
#SBATCH -t 02:00:00

# bbmerge Run Script
# This workflow sets a maximum length for the merged rbcL insert in order to output unmerged sequences.  
# These unmerged sequences resulting in too-long inserts will then filtered from the input to dada2. 

module load miniconda3
source activate /project/gbru_fy24_stinkbug_diet/annette/qiime_time # This should be deprecated but for some reason works better than conda activate

#Input variables

filtered="results_1453/filtered-seqs"
unmerged="results_1453/unmerged"
results="results_1453"
threads=3

# Data has already been trimmed and low-quality sequences removed.
# List all the forward reads 
file_list=$(find ${filtered} -type f -name "GT*L001_R1_001.fastq.gz")

# Trim primers from reads with bbtools
for file in $file_list; do
  dirname=$(dirname ${file})
  bname=$(basename ${file} L001_R1_001.fastq.gz)
  forward=${dirname}/${bname}L001_R1_001.fastq.gz
  plate=$(cut -d'_' -f1 <<<${bname})
  well=$(cut -d'_' -f2 <<<${bname})
  samp=$(cut -d'_' -f3 <<<${bname})
  tag=$(cut -d'_' -f4 <<<${bname})
  pattern=$plate"_"$well"_"$samp"_*_L001_R2_001.fastq.gz"
  reverse=$(find ${dirname}/ -name $pattern -printf "%p")  
  bbmerge.sh in=${forward} \
  in2=${reverse} \
  threads=${threads} \
  maxlength=135 \
  merge=f \
  outinsert=${results}/insert_sizes \
  ihist=${results}/insert_histogram \
  outu=${unmerged}/$(basename ${forward}) \
  outu2=${unmerged}/$(basename ${reverse}) \
  overwrite;
done
