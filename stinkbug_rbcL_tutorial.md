# Analysis of stinkbug diet chloroplast rbcL amplicon sequences

## Project outline

Chloroplast *rbcL* amplicon sequences were amplified from stinkbugs using the rbcLZ1 forward primer (ATGTCACCACAAACAGAGACTAAAGCAAGT) and the rbcL19 reverse primer (AGATTCCGCAGCCACTGCAGCCCCTGCTTC) (Poinar et al, 1998). This project processed 4701 samples from 8 Illumina Miseq runs: 
- Original run (2022)
- 1428
- 1450
- 1451
- 1452
- 1453
- 1454
- 1455

The sequences were demultiplexed either by hand (2022) or by the sequencing facility (all other runs).  The primers were trimmed, and sequences that had mismatched primer regions that could not be trimmed were filtered by length (135 bp).  Sequences that still had untrimmed primer regions were identified by testing the merged forward and reverse reads and filtering sequences again by length.  The clean sequences were then merged and amplicon sequence variants (ASV) were identified.  ASV tables from all eight sequencing runs were consolidated and taxonomically classified, and finally a robust Aitchison Principal Component Analysis (RPCA) was performed.  The RPCA showed that while there is much overlap among the sequencing runs, they do not fully overlap and further analysis might require correcting for batch effects.

## Demultiplexing
The 2022 samples were repaired using `repair.sh` from `BBTools` (BBMap – Bushnell B. – sourceforge.net/projects/bbmap/) and demultiplexed using `Ultraplex` (Wilkins, 2021) as described below and in `qiime_runscript.sh`. The rest of the data were demultiplexed by the sequencing facility. 

```
readdir="GT01-GT07_RBCL"
repaired="repaired"
demultiplexed="demultiplexed"
rbcl_f="ATGTCACCACAAACAGAGACTAAAGCAAGT"
rbcl_r="AGATTCCGCAGCCACTGCAGCCCCTGCTTC"
trimmed="fulltrimmed"
threads=3

# Combinatorial demultiplexing
# Build the barcode reference file
python data/ultraplex_barcodefiltermaker.py

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
```

## Trimming and filtering
After demultiplexing, primers were trimmed from the sequences using `bbduk.sh` from `BBTools` using a hamming distance of 2. Trimmed sequences were imported in `Qiime2` (2024.5; Bolyen et al, 2019) using a manifest and then summarized.  Note that `XXXX` denotes a four-digit ID for one of the sequencing runs. 

```
readdir="/project/gbru_fy24_stinkbug_diet/Illumina-XXXX"
rbcl_f="ATGTCACCACAAACAGAGACTAAAGCAAGT"
rbcl_r="AGATTCCGCAGCCACTGCAGCCCCTGCTTC"
trimmed="results_XXXX/fulltrimmed"
threads=3

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

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format 'PairedEndFastqManifestPhred33V2' \
  --input-path 'manifest_XXXX_filtered.tsv' \
  --output-path results_XXXX/input.qza

qiime demux summarize \
  --i-data results_XXXX/input.qza \
  --o-visualization results_XXXX/input.qzv
```

Some sequences weren't trimmed succesfully because of mismatches with the primers. The vast majority of sequences were 123 bp and `cutadapt` (Martin, 2011) was used to remove sequences that were longer than 135 bp.  The 2022 data had to be pre-processed to insert a space between the sample identifier and the "1" or "2" denoting forward and reverse reads before running `cutadapt`.

```
file_list=$(find results_2022/fulltrimmed -type f -name "ultraplex_GT*Fwd.fastq.gz") 
i=0
for file in $file_list; do
    i=$((i+1))
    if (( i % 100 == 0 )); then
        echo "File: $i"
    fi
    zcat $file | sed 's/1:N:0:/ 1:N:0:/g' | gzip > results_2022/fixed_read_names/$(basename $file)
done

file_list=$(find results_2022/fulltrimmed -type f -name "ultraplex_GT*Rev.fastq.gz") 
i=0
for file in $file_list; do
    i=$((i+1))
    if (( i % 100 == 0 )); then
        echo "File: $i"
    fi
    zcat $file | sed 's/2:N:0:/ 2:N:0:/g' | gzip > results_2022/fixed_read_names/$(basename $file)
done
```

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences results_XXXX/input.qza \
  --p-cores $threads \
  --p-max-n 135 \
  --o-trimmed-sequences results_XXXX/filtered-seqs.qza
```

Export the filtered reads `Qiime` artifact to fastq files.

```
qiime tools export \
    --input-path results_XXXX/filtered-seqs.qza \
    --output-path results_XXXX/filtered-seqs
```

We merged a sample pair of files from the trimmed data with `bbmerge.sh` (Bushnell, 2017) from `BBTools` and look at the average insert size to select Dada2 parameters

```
bbmerge.sh in=${trimmed}/GT01_S1_L001_R1_001.fastq.gz in2=${trimmed}/GT01_S1_L001_R2_001.fastq.gz
```

Some of the merged sequences were still too long due to untrimmed, mismatched primer sequences (160--210 bp).  A list of these sequences was made using `bbmerge.sh` by restricting the length of the sequences to 135 bp and saving only the unmerged sequences.

```
readdir="/project/gbru_fy24_stinkbug_diet/Illumina-XXXX"
filtered="results_XXXX/filtered-seqs"
unmerged="results_XXXX/unmerged"
results="results_XXXX"
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
```

The unmerged sequences were compared to the sequences filtered by `cutadapt` in Python, and those not included in the unmerged list were retained.

```
import os
import gzip
from Bio import SeqIO

def read_fastq_gz(file_path):
    """Read sequences from a gzipped FASTQ file and return a list of records."""
    records = []
    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            records.append(record)  # Store the entire record (including quality)
    return records

def find_unique_sequences(records1, records2):
    """Find unique sequences in two lists of FASTQ records based on IDs."""
    ids_file1 = {record.id: record for record in records1}
    ids_file2 = {record.id: record for record in records2}

    # Find unique sequences based on IDs
    unique_to_file1 = {id: ids_file1[id] for id in ids_file1 if id not in ids_file2}
    return unique_to_file1.values()
    
def save_sequences_to_fastq(records, output_file):
    """Save a list of FASTQ records to a FASTQ file."""
    with open(output_file, "w") as out_handle:
        SeqIO.write(records, out_handle, "fastq")  # Write in FASTQ format

def process_directories(dir1, dir2, output_dir):
    """Process matching FASTQ files in the specified directories and save unique sequences to the output directory."""
    fastq_files_dir1 = {f: os.path.join(dir1, f) for f in os.listdir(dir1) if f.endswith('.fastq.gz')}
    fastq_files_dir2 = {f: os.path.join(dir2, f) for f in os.listdir(dir2) if f.endswith('.fastq.gz')}
     
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Compare matching files in dir1 and dir2
    for file_name in fastq_files_dir1.keys():
        if file_name in fastq_files_dir2:
            path1 = fastq_files_dir1[file_name]
            path2 = fastq_files_dir2[file_name]
            print(f"Processing: {file_name}")

            records1 = read_fastq_gz(path1)
            records2 = read_fastq_gz(path2)

            unique_records1 = find_unique_sequences(records1, records2)

            # Save unique sequences to FASTQ files in the output directory
            save_sequences_to_fastq(unique_records1, os.path.join(output_dir, f"unique_seq_{file_name}"))

def main():
    # Specify the directories containing your fastq.gz files and the output directory
    dir1 = "results_XXXX/filtered-seqs"
    dir2 = "results_XXXX/unmerged"
    output_dir = "results_XXXX/to_dada"

    process_directories(dir1, dir2, output_dir)

if __name__ == "__main__":
    main()

```
The sequences to be sent to be used to construct the ASVs were imported as a `Qiime` artifact using a manifest file and summarized as a visualization.

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format 'PairedEndFastqManifestPhred33V2' \
  --input-path 'manifest_1452_to_dada.tsv' \
  --output-path results_1452/to_dada.qza

qiime demux summarize \
  --i-data results_1452/to_dada.qza \
  --o-visualization results_1452/to_dada.qzv
```
## Amplicon Sequence Variants

Run `dada2` (Callahan et al, 2016) to merge sequences and infer the ASVs. 
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs results_XXXX/to_dada.qza \
  --p-trunc-len-f 110 \
  --p-trunc-len-r 110 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-pooling-method pseudo \
  --p-n-threads $threads\
  --o-table results_XXXX/table_filtered.qza \
  --o-representative-sequences results_XXXX/repseqs_filtered.qza \
  --o-denoising-stats results_XXXX/dada2_stats_filtered.qza
```
Summarize the run.
```
qiime metadata tabulate \
  --m-input-file results_XXXX/dada2_stats_filtered.qza \
  --o-visualization results_XXXX/dada2_stats_filtered.qzv

```
Summarize the ASV table
```
qiime feature-table summarize \
  --i-table results_XXXX/table_filtered.qza \
  --o-visualization results_XXXX/table_filtered.qzv
```
## Taxonomy
### Reference Data

Our taxonomic classifier was based on the *rbcL* reference library DB4Q2 published by Dubois et al (2022). The library is available in [figshare](https://figshare.com/articles/online_resource/QIIME2_RefDB_development_zip/17040680?file=50431443).

```
qiime feature-classifier extract-reads \
    --i-sequences NCBI_rbcL_Viridiplantae_fasta_file_2021_06_14.qza \
    --p-f-primer $rbcl_f \
    --p-r-primer $rbcl_r \
    --p-n-jobs $threads \
    --o-reads rbcL_extracted.qza
```

There are 29,456 sequences in the extracted references.

### Naive Bayes Classifier

Train the classifier
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads rbcL_extracted.qza \
  --i-reference-taxonomy NCBI_rbcL_Viridiplantae_taxonomic_lineages_2021_06_14.qza \
  --o-classifier rbcL_classifier_extracted.qza
```

### Batch Analysis

Merge ASV tables together.
```
qiime feature-table merge \
  --i-tables results_2022/table_filtered.qza \
  --i-tables results_1428/table_filtered.qza \
  --i-tables results_1450/table_filtered.qza \
  --i-tables results_1451/table_filtered.qza \
  --i-tables results_1452/table_filtered.qza \
  --i-tables results_1453/table_filtered.qza \
  --i-tables results_1454/table_filtered.qza \
  --i-tables results_1455/table_filtered.qza \
  --o-merged-table merged_table_filtered.qza
```

Merge the representative sequence artifacts.
```
qiime feature-table merge-seqs \
  --i-data results_2022/repseqs_filtered.qza \
  --i-data results_1428/repseqs_filtered.qza \
  --i-data results_1450/repseqs_filtered.qza \
  --i-data results_1451/repseqs_filtered.qza \
  --i-data results_1452/repseqs_filtered.qza \
  --i-data results_1453/repseqs_filtered.qza \
  --i-data results_1454/repseqs_filtered.qza \
  --i-data results_1455/repseqs_filtered.qza \
  --o-merged-data merged_repseqs_filtered.qza
```

Merge metadata in Python.

```
import pandas as pd
import numpy as np

run_list=['2022', '1428', '1450', '1451', '1452', '1453', '1454', '1455']
big_meta = []

for seqrun in run_list:
    print(seqrun)

    file='manifest_' + seqrun + '_filtered.tsv'
    this_data=pd.read_csv(file, delimiter='\t')
    this_data['sequencing-run']=seqrun
    this_data = this_data.drop(columns = ['forward-absolute-filepath', 'reverse-absolute-filepath'])
    big_meta.append(this_data)
    
big_df = pd.concat(big_meta, ignore_index=True)
print(big_df.head())

big_df.to_csv("metadata_filtered_merged.tsv", sep = "\t", index=False, header=True)
```

Summarize merged artifact
```
qiime feature-table tabulate-seqs \
  --i-data merged_repseqs_filtered.qza \
  --o-visualization merged_repseqs_filtered.qzv

qiime feature-table summarize \
  --i-table merged_table_filtered.qza \
  --o-visualization merged_table_filtered.qzv \
  --m-sample-metadata-file metadata_filtered_merged.tsv
```

### Classify taxonomy 

Run the classifier on trimmed sequence data. 

```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier-trimmed.qza \
  --i-reads merged_repseqs_filtered.qza \
  --o-classification merged_taxonomy_output_filtered.qza
```

Create a visualization of the taxonomic classification.

```
qiime metadata tabulate \
  --m-input-file merged_taxonomy_output_filtered.qza \
  --o-visualization merged_taxonomy_output_filtered.qzv
```
Summarize the taxonomy data as barplots

```
qiime taxa barplot \
--i-table merged_table_filtered.qza \
--i-taxonomy merged_taxonomy_output_filtered.qza \
--m-metadata-file metadata_filtered_merged.tsv \
--o-visualization merged_taxa_bar_plots_filtered.qzv
```

## Principal component analysis

Perform Atchison's Robust PCA (RPCA) using package `gemelli` (Martino et al, 2019).  This package is not included in the `Qiime 2` amplicon suite, so I installed `gemelli` in a Conda environment.

```
conda activate qiime_time2; qiime gemelli rpca \
    --i-table merged_table_filtered.qza \
    --p-min-sample-count 500 \
    --o-biplot rbcl_ordination_filtered.qza \
    --o-distance-matrix rbcl_distance_filtered.qza
```

Build the emperor plots to visualize the RPCA (Vázquez-Baeza et al, 2013).
```
    qiime emperor biplot \
    --i-biplot rbcl_ordination_filtered.qza \
    --m-sample-metadata-file metadata_filtered_merged.tsv \
    --m-feature-metadata-file merged_taxonomy_output_filtered.qza  \
    --o-visualization rbcl_biplot_taxonomy_filtered.qzv
```

![Centered emporer plot](/images/RPCA_stinkbug_rbcl_centered.png)
![Rotated emporer plot](/images/RPCA_stinkbug_rbcl_rotated.png)

While the 8 runs mostly overlap each other, there appear to be two major clusters that make up each run to varying degrees.  For further analysis, we recommend controlling for batch effects.  


### Authors:  Annette M. Hynes and Adam R. Rivers

### References:


Bolyen, E., Rideout, J. R., Dillon, M. R., Bokulich, N. A., Abnet, C. C., Al-Ghalith, G. A., Alexander, H., Alm, E. J., Arumugam, M., Asnicar, F., Bai, Y., Bisanz, J. E., Bittinger, K., Brejnrod, A., Brislawn, C. J., Brown, C. T., Callahan, B. J., Caraballo-Rodríguez, A. M., Chase, J., … Caporaso, J. G. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, **37**(8), 852–857. https://doi.org/10.1038/s41587-019-0209-9

B. Bushnell, J. Rood, and E. Singer. 2017. BBMerge – Accurate paired shotgun read merging via overlap . *PLoS ONE*. **12**(10) : e0185056. 

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nature Methods*, **13**(7), 581. https://doi.org/10.1038/nmeth.3869

Coordinators, N. R. (2017). Database Resources of the National Center for Biotechnology Information. *Nucleic Acids Research*, **45**(D1), D12–D17. https://doi.org/10.1093/nar/gkw1071

Dubois et al (2022, *Genomic Data* **23**:53 doi: 10.1186/s12863-022-01067-5)

Johnson, M., Zaretskaya, I., Raytselis, Y., Merezhuk, Y., McGinnis, S., & Madden, T. L. (2008). NCBI BLAST: a better web interface. *Nucleic Acids Research*, **36**(suppl_2), W5–W9. https://doi.org/10.1093/nar/gkn201

Martino, C. et al (2019). A Novel Sparse Compositional Technique Reveals Microbial Perturbations. *mSystems* **4**, 

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet. Journal*, **17**(1), pp-10. https://doi.org/10.14806/ej.17.1.200

McKinney W. (2010). Data Structures for Statistical Computing in Python. In S. van der Walt & Jarrod Millman (Eds.), *Proceedings of the 9th Python in Science Conference* (pp. 51–56).

Poinar H.N., Hofreiter M., Spaulding W.G., Martin P.S., Stankiewicz B.A., Bland H., Evershed R.P., Possnert G., Pääbo S. (1998). Molecular Coproscopy: dung and diet of the extinct ground Sloth *Nothrotheriops shastensis*. *Science*. **281**:402–406. doi: 10.1126/science.281.5375.402

Vázquez-Baeza, Y., Pirrung, M., Gonzalez, A., & Knight, R. (2013). EMPeror: a tool for visualizing high-throughput microbial community data. *Gigascience*, **2**(1), 16. https://doi.org/10.1186/2047-217X-2-16

Wilkins O.G., Capitanchik C., Luscombe N.M., Jernej Ule. (2021). Ultraplex: A rapid, flexible, all-in-one fastq demultiplexer. *Wellcome Open Res*, **6**:141. doi: 10.12688/wellcomeopenres.16791.1