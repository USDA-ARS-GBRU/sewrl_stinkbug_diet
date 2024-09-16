# Analysis of Amplicon sequences of the RbcL chloroplast gene Locus

Adam Rivers  
October 17, 2023  
USDA-ARS-GBRU

# Project outline
This project processed 672 samples. The samples were multiplexed with combinatorial inline barcodes   
and then second indexed with Illumina barcodes. The data were sequences in 2 x 151 mode on an Illumina Miseq

# Data processing
The data were processed with bbtools, Ultraplex and Qiime2.
the script `qiime_runscript.sh` provides a general overview of the workflow.  

# Classifier Training

The RbcL reference dataset used for Qiime2 classifier training is from  https://doi.org/10.6084/m9.figshare.c.5504193.v1

Bell KL, Loeffler VM, Brosi BJ. An rbcL reference library to aid in the identification 
of plant species mixtures by DNA metabarcoding. Appl Plant Sci. 2017 Mar 10;5(3):apps.1600110.
doi: 10.3732/apps.1600110. 

* This script splits the data in the Bell dataset into a two-column taxonomy.txt file and a fasta file.
* Sequences are identified by an md5 hash rather than a taxonomy string with gaps that make the ids unreliable
* Taxids are removed from the taxonomy string for use in Qiime2
* Reformats some randomly interspersed RNA sequences to DNA
* Removes exact sequence duplicates in from the reference data set, retaining the first entry encountered. An RNA 
    sequence transcribed from the same DNA sequence is also considered identical and removed.

this was processed with the script `rbcl_db_parser.py`

