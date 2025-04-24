# Analysis of Amplicon sequences of the RbcL chloroplast gene Locus

Annette Hynes and Adam Rivers  
April 24, 2025  
USDA-ARS-GBRU

## Project outline
This project processed 4701 samples from 8 runs. The samples were multiplexed with combinatorial inline barcodes   
and then second indexed with Illumina barcodes. The data were sequences in 2 x 151 mode on an Illumina Miseq

## Data processing
The data were processed with bbtools, Ultraplex and Qiime2. The document `stinkbug_rbcL_tutorial.md` provides a general overview of the workflow and provides the code used.  

## Classifier training
The RbcL reference dataset used for Qiime2 classifier training is from Dubois et al (2022). The library is available in [figshare](https://figshare.com/articles/online_resource/QIIME2_RefDB_development_zip/17040680?file=50431443).

## Output files

## References 
Dubois, B., Debode, F., Hautier, L. et al.(2022). A detailed workflow to develop QIIME2-formatted reference databases for taxonomic analysis of DNA metabarcoding data. *BMC Genom Data* **23**:53  https://doi.org/10.1186/s12863-022-01067-5


