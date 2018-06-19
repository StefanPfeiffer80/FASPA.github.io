# FASPA - Fast Amplicon Sequence Processing and Analysis for MiSeq paired end sequence data

If you have any questions, please write me at pfeiffer.stefan@gmx.at
FASP is a workflow for analysing Illumina paired-end sequence data. 
This file is distributed without warranty
This file runs as a Linux bashscript; 
Cite as: Pfeiffer, S. (2018) Prootocl Exchange. FASPA - Fast Amplicon Sequence Processing and Analyses

# Introduction
High-throughput sequencing of the 16S rRNA marker gene is the current benchmark in the characterization of bacterial microbial communities from virtually all environments.
Here, I present FASPA (Fast Amplicon Sequence Processing and Analysis), an amplicon processing workflow that is easy to use and allows in depth-analysis of microbial communities. 
FASPA addresses the need for a transparent pipeline that based on executable bash scripts which allows very fast processing (less than 30 minutes on an average speed laptop with 4 GB RAM).
FASPA utilizes well-known and state of the art bioinformatics tools and provieds full transparency of the tools used and how they are applied.
FASPA also supports the integration of the processed amplicon data in various popular analysis tools, such as QIIME or the R-based (cite) package collection phyloseq via scripts that 
These pipelines claim that their usage addresses the need to be easy-to-use through the application of default or streamlined parameters. 
In a nutshell, FASPA manages precarious balance by being very fast, applies state of the art bioinformatic tools, having low CPU requirements, offers full transparency for the user by being at the same time easy to use. 

# Main structure
# Executable UNIX shell bash scripts
# *FASP_preprocess.sh, FASP_unoise.sh, FASP_uparse.sh*
Put the bash scripts in the folder where your fastq files are. Open the script in a text editor. IMPORTANT!!!!!!! Positions marked with XX need to be adjusted according to the users need!!!!!!!
1. FASP_Preprocessing.sh
Bash script for the preprocessing of raw fastq files based on the programs USEARCH v10.240 and VSEARCH v2.80 !links - check markdown!!!.
What the script does: 
a. Merging forward and reverse paired-end reads into one single .fastq file -> more information http://drive5.com/usearch/manual/merge_pair.html. Further, an info file on the single merged file is created, named raw_info.txt
b. Estimating the error rates bases based on of maybe screenshot?
![GitHub Logo](/logo.png)
Expected errors:
http://drive5.com/usearch/manual/exp_errs.html -> Output file is: *qualrawfq.txt*
c. Trimming of primers, overhangs and quality filtering of the reads based on estimated error rates.
d. Generation of a fasta file containing only unique sequences


Transparency: FASPA is a completely transparent workflow, advantegous for the user, also gives respect to the used programs.
4.	SINTAX algorithm for the taxonomic assignment and classification of OTUs or ASVs.
5.	Custom perl and R scripts that allow downstream analysis of the generated dataset i downstream analysis to utilize the statistical features of the popular tools QIIMEX the  R-based package phyloseqX and the Rhea-pipeline (L) which allows 


# 
