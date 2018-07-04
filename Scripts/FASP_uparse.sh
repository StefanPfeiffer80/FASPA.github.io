# FASP_uparse file; written by Stefan Pfeiffer, 1.December 2017, last modified on 30.6.2018; 
# Contact: microbiawesome@gmail.com; 
# Copyright (C) Stefan Pfeiffer, 2016-2018, all rights reserved.
# FASP is a workflow for analysing Illumina paired-end sequence data. 
# IMPORTANT!!!!!!! You can adjust positions manually in the script and safe the script under a new file name.
# This file is distributed without warranty
# This file runs as a Linux bashscript; You can run this script only after finishing FASP_preprocessing.sh;
# This script generates an OTU table, a phylogenetic tree, a taxonomy file
# Cite as: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses, DOI:10.5281/zenodo.1302799

# 1.Cluster OTUs with UPARSE, at 97% sequence similarity 
./US_10.240 -cluster_otus uniques.fa -otus otus.fa -relabel Otu
echo "!!!!!"
echo "OTU clustering successfull; Output file: otus.fa"
echo "!!!!!"
# 2. Sequence length trimming
echo "Sequence length trimming, did you check the expected length of your OTUs?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done
while getopts ":i:" opt
   do
     case $opt in
        i ) min=$OPTARG;;

     esac
done
./US_10.240 -sortbylength otus.fa -fastaout otus_sorted_UP.fasta -minseqlength $min # Cut all OTUs that are shorter than a minimum sequence length!!!
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Sequence length trimming successfull; Output file: otus_sorted_UP.fasta"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#3. Make TU table
./US_10.240 -otutab raw.fq -otus otus_sorted_UP.fasta -otutabout otutab_UP_raw.txt
#4. UncrossOTU table to remove spurious OTUs generated via cross-talk 
./US_10.240 -uncross otutab_UP_raw.txt -tabbedout outUP.txt -report repUP.txt -otutabout otutab_UP_uncrossed.txt
#6. Statistics of your OTU table
./US_10.240 -otutab_stats otutab_UP_uncrossed.txt -output otutab_UP_uncrossed_STATS.txt 
#7. Assign taxonomy using SINTAX
# Default value for the sintax_cutoff is 0.5, can be adjusted
./US_10.240 -sintax otus_sorted_UP.fasta -db rdp_16s_v16.fa -strand both \-tabbedout sintaxotusrdp.txt -sintax_cutoff 0.5 
#8. Delete OTUs from the SINTAX-OUTPUT FILE file that are not in the OTU table
perl ./project/genomics/Stefan/FASPA/FASP_tax_filtered.pl otutab_UP_uncrossed.txt sintaxotusrdp.txt SINTAX_OTUS_RDP_FILT.txt
#9. Make a phylogenetic tree
./US_10.240 -cluster_agg otus_sorted_UP.fasta -treeout Tree_UP.tree
#10.Construct a QIIME OTU table .txt and BIOM using Greengenes annotation -> if you want to use PICRUSt for example
perl ./Otu_tab_tax_greengenes.pl otutab_UP_uncrossed.txt SINTAX_OTUS_RDP_FILT.txt otutab_otus_UP_rdp_greengenes_format.txt
#10.Construct a QIIME OTU table .txt and BIOM using SILVA annotation
perl ./Otu_tab_tax_SILVA.pl otutab_UP_uncrossed.txt SINTAX_OTUS_RDP_FILT.txt otutab_otus_UP_rdp_SILVA_format.txt

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!!!Amplicon proscessing is completed!!!!!!!!!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "For the R-workflow you will need the following files:"
echo "SINTAX_OTUS_RDP_FILT  -> R - TAX-file"
echo "otutab_UP_uncrossed   -> R - OTU-file"
echo "Tree_UP.tree          -> R - Tree-file"
