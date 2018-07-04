# FASP peprocessing file; written by Stefan Pfeiffer, 1.December 2017, last modified on 4.7.2018; 
# Contact: microbiawesome@gmail.com; 
# Copyright (C) Stefan Pfeiffer, 2016-2018, all rights reserved.
# FASPA is a workflow for analysing Illumina paired-end sequence data. 
# This file is distributed without warranty
# This file runs as a Linux bashscript; 
# Cite as: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses, DOI:10.5281/zenodo.1302800

set -e
set -o pipefail

# 1.Merge reads: All read pairs will be merged into a single file
./US_10.240 -fastq_mergepairs *_R1_001.fastq \-fastqout raw.fq -relabel @ 
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Read pairs are merged successfully!!!! Output file: raw.fq"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 2. Raw_quality_info_merged_reads: Gives an output on the filesize, the quality and the sequence length -> Output file: "raw_info.txt"
./US_10.240  -fastx_info raw.fq -secs 5 -output raw_info.txt  
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Successful; You get information on the filesize of raw.fq!!!!!"
echo "length of sequences, base frequencies and read qualitities!!!!"
echo "!Output file: raw_info.fq!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 4. Quality check of the sequences based on eestats (expected error calculation)
./US_10.240  -fastq_eestats2 raw.fq -output qualrawfq.txt #Output file: "qualrawfq.txt"; a table that provides you the maximum expected error rates for a sequence - e.g. MaxEE = expected error is 0 -> see Edgar & Flyvberg 2014
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Successful; Output file: raw_info.fq"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 5. Find out at which positions your primers are
./US_10.240  -fastx_subsample raw.fq -sample_size 100 -fastqout raw_subset_100.fq # The subset sample size can be changed, default is 100
./US_10.240  -search_oligodb raw_subset_100.fq -db primers.fa -strand both -userout primer_positions.txt -userfields query+qlo+qhi+qstrand
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Did you check the length of your primers and the expected size of your amplicons?"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done
# 6. Removal of primers + adapters
while getopts ":l:r:" opt
   do
     case $opt in
        l ) primerleft=$OPTARG;;
        r ) primerright=$OPTARG;;
        q ) ee=$OPTARG;;

     esac
done
./US_10.240 -fastx_truncate raw.fq -stripleft $primerleft -stripright $primerright \-fastqout strippedraw.fq # Output file: "strippedraw.fq";  left=forward primer; right = reverse primer; e.g. -stripleft 20 -stripright 19
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Successful, Primers and Adaptors are removed; Output file: strippedraw.fq"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 7.Quality filtering of sequences: We use the number of expected errors (=sum of error probabilities) as the benchmark of quality filtering
./US_10.240  -fastq_filter strippedraw.fq -fastq_maxee $ee \-fastaout filteredstripped.fa -relabel Filt # Output file: "filteredstripped.fa"; expected error rates can be set -fastq_maxee e.g. 1.0
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Successful, low quality reads were removed; Output file: filteredstripped.fa"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# 8.Extracting uniques sequences: 
./US_10.240  -fastx_uniques filteredstripped.fa -sizeout -relabel Uniq -fastaout uniques.fa # output file is "uniques.fa"
echo "Successful; The file uniques.fa contains your unique read seuqences; Next comes denoising (bash ./FASP_unoise.sh)or OTU clusetring (bash FASP_cluster.sh)"
