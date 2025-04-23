#!/bin/bash 

#SBATCH --mail-type=END,FAIL
#SBATCH --account=gbru_fy24_stinkbug_diet
#SBATCH --job-name=bbduk_1428
#SBATCH --output=bbduk_1428_%j.out

# bbduk Run Script
# This workflow trims primer sequences from Amplicon Sequence Variants 
# of RbcL amplicons with bbduk. I removed the error-correcting step because it
# caused downstream PHRED quality value errors. 

module load miniconda3
source activate /project/gbru_fy24_stinkbug_diet/annette/qiime_time # This should be deprecated but for some reason works better than conda activate

#Input variables

readdir="/project/gbru_fy24_stinkbug_diet/Illumina-1428"
rbcl_f="ATGTCACCACAAACAGAGACTAAAGCAAGT"
rbcl_r="AGATTCCGCAGCCACTGCAGCCCCTGCTTC"
trimmed="results_1428/fulltrimmed"
threads=3

# Data has already been demultiplexed
# List all the forward reads without any "UNKNOWN" files 
file_list=$(find ${readdir} -type f -name "LJLN8_s1_R1*_GT*.fastq.gz")

# Trim primers from reads with bbtools
for file in $file_list; do
  dirname=$(dirname ${file})
  bname=$(basename ${file})
  prename=${bname#'LJLN8_s1_R1'}
  forward=${dirname}/LJLN8_s1_R1${prename}
  reverse=${dirname}/LJLN8_s1_R2${prename}
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

