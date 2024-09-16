#!/bin/bash 
# Qiime Run Script
# This workflow calls Amplicon Sequence Variants of RbcL amplicons with Qiime2 and Dada2
# and trains a Niaee Bayes taxonomic classifier  

#activate qiime2 amplicon Conda environment 

# conda activate qiime2-amplicon-2023.9

#Input varaibles

readdir="GT01-GT07_RBCL"
repaired="repaired"
demultiplexed="demultiplexed"
rbcl_f="ATGTCACCACAAACAGAGACTAAAGCAAGT"
rbcl_r="AGATTCCGCAGCCACTGCAGCCCCTGCTTC"
trimmed="fulltrimmed"
threads=3

# Combinatorial demultiplexing
# python data/ultraplex_barcodefiltermaker.py

# repair disordered reads 

for file in ${readdir}/*R1_001.fastq.gz; do
    bname=`basename $file _R1_001.fastq.gz`
    forward=${bname}_R1_001.fastq.gz
    reverse=${bname}_R2_001.fastq.gz
    repair.sh in=${readdir}/${forward} \
    in2=${readdir}/${reverse} \
    out=${repaired}/${forward} \
    out2=${repaired}/${reverse} \
    tossbrokenreads ;
  done

# demultiplex files
for file in ${repaired}/*R1_001.fastq.gz; do
    bname=`basename $file _R1_001.fastq.gz`
    forward=${bname}_R1_001.fastq.gz
    reverse=${bname}_R2_001.fastq.gz
    ultraplex \
    -i ${repaired}/${forward} \
    -i2 ${repaired}/${reverse} \
    -b data/ultraplex_barcodes.csv \
    --dont_build_reference \
    -m3 1 \
    -inm \
    -t $threads \
    -d ${demultiplexed} \
    -o $bname
done




# First trim barcodes from the reads with bbtools
for file in ${demultiplexed}/*_Fwd.fastq.gz; do
    bname=`basename $file _Fwd.fastq.gz`
    forward=${bname}_Fwd.fastq.gz
    reverse=${bname}_Rev.fastq.gz
    bbduk.sh in=${demultiplexed}/$forward \
         in2=${demultiplexed}/${reverse} \
         tossbrokenreads \
         barcodes \
         threads=${threads} \
         literal=${rbcl_f},${rbcl_r} \
         out=${trimmed}/${forward} \
         out2=${trimmed}/${reverse} \
         overwrite \
         ktrim=rl;
done

# Import reads as Qiime artifact
# names must be in format 'GT02_S2_L001_R2_001.fastq.gz'
#qiime tools import \
#  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-path $trimmed \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --output-path input.qza

# python snippet to find some samples that were missing
python manifest_filter.py > data/manifest_fliltered.txt


qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data/manifest_filtered.txt \
  --output-path input.qza \
  --input-format PairedEndFastqManifestPhred33V2



# Visualize trimmed reads
qiime demux summarize --i-data input.qza --o-visualization input.qzv
# qual scores are high up to about 110

# We will merge trimmed data with BBerge and look at the average insert size dto select Dada2 parameters

#bbmerge.sh in=trimmed/GT01_S1_L001_R1_001.fastq.gz in2=trimmed/GT01_S1_L001_R2_001.fastq.gz

# Avg Insert:          	124.6
# Standard Deviation:  	9.2
# Mode:                	123

# Insert range:        	107 - 204
# 90th percentile:     	123
# 75th percentile:     	123
# 50th percentile:     	123
# 25th percentile:     	123
# 10th percentile:     	118
# even at 200 nt a 130+130 fragment would give 60 bp of overap, more than the needed 20 bp + biolgocial variation 

# 
# Run Data2 to create the ASVs 
qiime dada2 denoise-paired --i-demultiplexed-seqs input.qza \
                    --p-trunc-len-f 110 \
                    --p-trunc-len-r 110 \
                    --p-trim-left-f 0 \
                    --p-trim-left-r 0 \
                    --p-pooling-method pseudo \
                    --p-n-threads $threads\
                    --o-table table.qza \
                    --o-representative-sequences repseqs.qza \
                    --o-denoising-stats dada2_stats.qza

#Summarize the run
qiime metadata tabulate --m-input-file dada2_stats.qza --o-visualization dada2_stats.qzv

# Summarize the ASV table
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv


# Prepare Reference data
# We use the RbcL reference from: 
# Bell KL, Loeffler VM, Brosi BJ. An rbcL reference library to aid in the identification 
# of plant species mixtures by DNA metabarcoding. Appl Plant Sci. 2017 Mar 10;5(3):apps.1600110.
# doi: 10.3732/apps.1600110. 
# Download RbcL reference dataset from  https://doi.org/10.6084/m9.figshare.c.5504193.v1extract and place in "reference/"
# This data seet needs to be cleaned and parsed  with my script rbcl_db_parser.py prior to use

# python rbcl_db_parser.py

# #import The cleaned RbcL reference fasta file Qiime2
# qiime tools import --input-path renamed.fa \
#     --output-path refseqs.qza \
#     --type 'FeatureData[Sequence]'

# # import the new Taxonomy file into Qiime2
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-format HeaderlessTSVTaxonomyFormat \
#   --input-path taxonomy.txt \
#   --output-path ref-taxonomy.qza

# # Trim reads 
# qiime feature-classifier extract-reads \
#   --i-sequences refseqs.qza \
#   --p-f-primer $rbcl_f \
#   --p-r-primer $rbcl_r \
#   --p-min-length 100 \
#   --p-max-length 200 \
#   --p-n-jobs $threads \
#   --o-reads refseqs_trimmed.qza


# # Train the Naive Bayes Classifier

# qiime feature-classifier fit-classifier-naive-bayes \
#   --i-reference-reads refseqs_trimmed.qza \
#   --i-reference-taxonomy ref-taxonomy.qza \
#   --o-classifier classifier-trimmed.qza

# Run the classifier on trimmed classifier data 

  qiime feature-classifier classify-sklearn \
  --i-classifier classifier-trimmed.qza \
  --i-reads repseqs.qza \
  --o-classification taxonomy_outout-trimmed.qza

# Create a visualization of the taxonomic classification

qiime metadata tabulate \
  --m-input-file taxonomy_outout-trimmed.qza \
  --o-visualization taxonomy_outout-trimmed.qzv

# summarize the data as barplots

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy_outout-trimmed.qza \
--m-metadata-file data/manifest.txt \
--o-visualization taxa-bar-plots.qzv

