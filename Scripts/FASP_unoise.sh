# FASP peprocessing file; written by Stefan Pfeiffer, 1.December 2017,  last modified on 30.6.2018. 
# Contact: microbiawesome@gmail.com; 
# Copyright (C) Stefan Pfeiffer, 2016-2018, all rights reserved.
# This file is distributed without warranty
# This file runs as a Linux bashscript; You can run this script only after finishing FASP_preprocessing.sh;
# This script generates an OTU table, a phylogenetic tree, a taxonomy file
# Cite as: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses, DOI:10.5281/yenodo.1302800

# 1.Denoising of reads -> 
./US_10.240 -unoise3 uniques.fa -zotus zotus.fa
# 2. Change ZOTUs(Zero-radius OTUs) to OTUs in the fasta file to enable downstream processing
sed -i "s/\Zo/\O/g" zotus.fa
echo "!!!!!"
echo "Denoising successfull; Output file: Zotus.fa"
echo "!!!!!"
# 3. Sequence length trimming
echo "Sequence length trimming, did you check the expected length of your ZOTUs/OTUs?"
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
./US_10.240 -sortbylength zotus.fa -fastaout otus_sorted_UN.fasta -minseqlength $min # Cut all sequences that are shorter than a minimum sequence length!!!
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Sequence length trimming successfull; Output file: otus_sorted_UN.fasta"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#4. Make (Z)OTU table
./US_10.240 -otutab raw.fq -otus otus_sorted_UN.fasta -otutabout otutab_UN_raw.txt
#5. Uncross(Z)OTU table to remove spurious (Z)OTUs generated via cross-talk 
./US_10.240 -uncross otutab_UN_raw.txt -tabbedout outUN.txt -report repUN.txt -otutabout otutab_UN_uncrossed.txt
#6. Statistics of your (Z)OTU table
./US_10.240 -otutab_stats otutab_UN_uncrossed.txt -output otutab_UN_uncrossed_STATS.txt 
#7. Assign taxonomy using SINTAX
./US_10.240 -sintax otus_sorted_UN.fasta -db rdp_16s_v16.fa -strand both \-tabbedout sintaxzotusrdp.txt -sintax_cutoff 0.5 # Default is 0.5
#8. Delete (Z)OTUs from the SINTAX-OUTPUT FILE file that are not in the OTU table
perl ./FASP_tax_filtered.pl otutab_UN_uncrossed.txt sintaxzotusrdp.txt SINTAX_ZOTUS_RDP_FILT.txt
#9. Make a phylogenetic tree
./US_10.240 -cluster_agg otus_sorted_UN.fasta -treeout Tree_UN.tree
#10.Construct a QIIME OTU table .txt and BIOM using Greengenes annotation -> if you want to use PICRUSt for example
perl ./Otu_tab_tax_greengenes.pl otutab_UN_uncrossed.txt SINTAX_ZOTUS_RDP_FILT.txt otutab_otus_UN_rdp_greengenes_format.txt
#10.Construct a QIIME OTU table .txt and BIOM using SILVA annotation
perl ./Otu_tab_tax_SILVA.pl otutab_UN_uncrossed.txt SINTAX_ZOTUS_RDP_FILT.txt otutab_otus_UN_rdp_SILVA_format.txt

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "!!!!!!!!!Amplicon proscessing is completed!!!!!!!!!!!"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "For the R-workflow you will need the following files:"
echo "SINTAX_ZOTUS_RDP_FILT -> R - TAX-file"
echo "otutab_UN_uncrossed   -> R - OTU-file"
echo "Tree_UN.tree             R - Tree-file"
