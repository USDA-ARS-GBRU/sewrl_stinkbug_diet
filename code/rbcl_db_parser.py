"""rbcl_db_parser: a script to reformat and correct errors in the 
RbcL reference dataset from  https://doi.org/10.6084/m9.figshare.c.5504193.v1for use in Qiime2.

Adam Rivers. USDA-ARS. October 13, 2023 

Bell KL, Loeffler VM, Brosi BJ. An rbcL reference library to aid in the identification 
of plant species mixtures by DNA metabarcoding. Appl Plant Sci. 2017 Mar 10;5(3):apps.1600110.
doi: 10.3732/apps.1600110. 

Changes the  Script makes to the dataset:

* This script splits the data in the Bell dataset into a two column taxonomy.txt file and a fasta file.
* Sequences are idenitifed by an md5 hash rather than a taxonomy string with gaps that makes the ids unreliable
* Taxids are removed from the taxonomy string for use in Qiime2
* Reformats some randomly interspersed RNA sequences to DNA
* Removes exact sequence duplicates in from the reference data set, retaining the first entry encountered. An RNA 
    sequence transcribed from the same DNA sequence are also considered identical and removed.


"""

import hashlib

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas

#todo check for RNA
# input
#>k__Viridiplantae_33090;p__Streptophyta_35493;c__sub__asterids_71274;o__Solanales_4069;f__Solanaceae_4070;g__Grabowskia_84167;s__Grabowskia glauca_84168;
# level__levename__taxid;
#needed output
# 229854	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
#>1111561 AGAGTTTGATCCTGGCTCAGAT
reffile = "reference/rbcL_plus_Nov2019adds_Jul2021corrections.dada2.fa"
#reffile = "reference/sample.fa" 
dada2in = SeqIO.parse(reffile, 'fasta') 
seqid = []
taxonomy = []
hashlist = []
with open("renamed.fa", 'w') as fastaout:
    for record in dada2in:
        record.seq = record.seq.back_transcribe()
        md5 = hashlib.md5(str(record.seq).encode()).hexdigest()
        if md5 not in hashlist:
            hashlist.append(md5)
            long = record.description
            longlist = long.split(';')
            cleanlist = []
            for item in longlist:
                a = item.split("_")
                cleanlist.append("_".join(a[:-1]))
            outline = ";".join(cleanlist)
            md5 = hashlib.md5(str(record.seq).encode()).hexdigest()
            seqid.append(md5)
            taxonomy.append(outline)
            record.id = md5
            SeqIO.write(record, fastaout, "fasta")


taxdf = pandas.DataFrame.from_dict({'id': seqid, 'taxonomy': taxonomy})
taxdf.to_csv("taxonomy.txt", header=None, index=None, sep='\t')