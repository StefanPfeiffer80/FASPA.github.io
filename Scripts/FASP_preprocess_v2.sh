# Contact: microbiawesome@gmail.com; 
# Copyright (C) Stefan Pfeiffer, 2016-2018, all rights reserved.
# FASP is a workflow for analysing Illumina paired-end sequence data. 
# This file is distributed without warranty
# This file runs as a Linux bashscript; 
# Cite as: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses, DOI:10.5281/zenodo.1302799

set -e
set -o pipefail

# 1.Merge reads: All read pairs will be merged into a single file
./US_10.240 -fastq_mergepairs *_R1_001.fastq \-fastqout raw.fq -relabel @ 
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Read pairs are merged successfully!!!! Output file: raw.fq"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 2. Find out at which positions your primers are
vsearch  -fastx_subsample raw.fq -sample_size 100 -fastqout raw_subset_100.fq # The subset sample size can be changed, default is 100
./US_10.240  -search_oligodb raw_subset_100.fq -db primers.fa -strand both -userout primer_positions_v2.txt -userfields query+qlo+qhi+qstrand
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Did you check the length of your primers and the expected size of your amplicons?"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done
# 3. Removal of primers + adapters
while getopts ":l:r:m:s:q:" opt
   do
     case $opt in
        l ) primerleft=$OPTARG;;
        r ) primerright=$OPTARG;;
	m ) max=$OPTARG;;
	s ) min=$OPTARG;;
  q ) ee=$OPTARG;;
     esac
done
### You have to enter minimum length of your sequences and the maximum lenght of your sequences. Due to differences in 16S amplicon length you should use a frame of 50-100 base pair positions (e.g. 300-360 for an amplicon of the expected size of 330 bp). 
vsearch -fastq_filter raw.fq --fastq_stripleft $primerleft --fastq_stripright $primerright --fastq_maxee $ee --fastq_maxlen $max --fastq_minlen $min --fastaout fileredstripped.fa
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Successful, low quality reads were removed; Output file: filteredstripped.fa"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 4.Extracting uniques sequences: 
./US_10.240 -fastx_uniques filteredstripped.fa -sizeout -relabel Uniq -fastaout uniques.fa # output file is "uniques.fa"
echo "Successful; The file uniques.fa contains your unique read seuqences; Next comes denoising (bash FASP_unoise.sh)or OTU clusetring (bash FASP_uparse.sh)"
