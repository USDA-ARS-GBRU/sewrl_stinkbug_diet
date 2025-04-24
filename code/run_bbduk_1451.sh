#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --account=gbru_fy24_stinkbug_diet
#SBATCH --job-name=bbduk_1451
#SBATCH --output=bbduk_1451_%j.out

# bbduk Run Script
# This workflow trims primer sequences from Amplicon Sequence Variants 
# of RbcL amplicons with bbduk. I removed the error-correcting step because it
# caused downstream PHRED quality value errors. 

module load miniconda3
source activate /project/gbru_fy24_stinkbug_diet/annette/qiime_time # This should be deprecated but for some reason works better than conda activate

#Input variables

readdir="/project/gbru_fy24_stinkbug_diet/Illumina-1451"
rbcl_f="ATGTCACCACAAACAGAGACTAAAGCAAGT"
rbcl_r="AGATTCCGCAGCCACTGCAGCCCCTGCTTC"
trimmed="results_1451/fulltrimmed"
threads=3

# Data has already been demultiplexed
# List all the forward reads without any "UNKNOWN" files 
file_list=$(find ${readdir} -type f -name "GT*.R1.fastq.gz")

# Trim primers from reads with bbtools
for file in $file_list; do
  dirname=$(dirname ${file})
  bname=$(basename ${file} .R1.fastq.gz)
  forward=${dirname}/${bname}.R1.fastq.gz
  reverse=${dirname}/${bname}.R2.fastq.gz
  bbduk.sh in=${forward} \
  in2=${reverse} \
  tossbrokenreads \
  threads=${threads} \
  literal=${rbcl_f},${rbcl_r} \
  ktrim=rl \
  k=23 \
  mink=11 \
  hdist=2 \
  tpe=t \
  tbo=t \
  maq=10 \
  out=${trimmed}/$(basename ${forward}) \
  out2=${trimmed}/$(basename ${reverse}) \
  overwrite;
done
