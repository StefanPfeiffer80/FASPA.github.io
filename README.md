

# FASPA - Fast Amplicon Sequence Processing and Analysis for MiSeq paired end sequence data

If you have any questions, please write me at pfeiffer.stefan@gmx.at.  
FASPA is a workflow for analysing Illumina paired-end sequence data. 
FASPA is a collection of shell bash scripts, perl scripts and R scripts that interact with state of the art programs used in sequence processing, USEARCH and VSEARCH
FASPA is distributed without warranty.  
Cite as: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses. DOI:kommt noch  

# Introduction
High-throughput sequencing of the 16S rRNA marker gene is the current benchmark in the characterization of bacterial microbial communities from virtually all environments.  
Here, I present FASPA -Fast Amplicon Sequence Processing and Analysis, an amplicon processing workflow that is easy to use and allows in depth-analysis of microbial communities. 
FASPA addresses the need for a transparent pipeline that based on executable bash scripts, perl scripts and  which allows very fast processing (less than 30 minutes on an average speed laptop with 4 GB RAM).
FASPA utilizes well-known and state of the art bioinformatics tools and provieds full transparency of the tools used and how they are applied.
FASPA also supports the integration of the processed amplicon data in various popular analysis tools, such as QIIME or the R-based (cite) package collection phyloseq via scripts that 
These pipelines claim that their usage addresses the need to be easy-to-use through the application of default or streamlined parameters. 
In a nutshell, FASPA manages precarious balance by being very fast, applies state of the art bioinformatic tools, having low CPU requirements, offers full transparency for the user by being at the same time easy to use. 

# Software needed
For the first part of the FASPA workflow, the amplicon processing, FASPA uses the programs USEARCH v.10.240 and optionally VSEARCH v.2.80. 
Usearch by Robert Edgar (Edgar 2010) is a collection of functions and algorithms to efficiently, fast and accurately transform raw amplicon reads into an OTU table for downstream analysis. Usearch can be downloaded as a single executable file (www.drive5.com/usearch/download.html). Usearch includes the uparse algorithm (Edgar 2013) to cluster OTUs, which showed improved accuracy in OTU assignment towards many other commonly used clustering algorithms and was already cited several thousand times. For details how UPARSE works, see here (LINK). Recently, Usearch implemented unoise (unoise paper), an algorithm for denoising of raw sequences, which actually means that the genetic variation of sequences is analyzed to find out what causes the sequence variation: Real sequence differences or sequencing errors. Nowadays, denoising of raw amplicon reads becomes more popular and is recommended by several experts in the field (citations). Usearch however, is free of charge only in the 32 bit version, which holds a 4GB memory cap. While this is not a problem for most datasets (depending on the sample type between 50 and 100 samples can be processed with the 32bit version), larger datasets will not be working out. For this reason, FASPA includes Vsearch by Torbjørn Rognes (Rognes et al. 2016), which was designed as an open source alternative to usearch. Vsearch can be downloaded here https://github.com/torognes/vsearch. For the FASPA workflow, it is evident that you don't use any vsearch version prior to v2.8.0.

**Usearch installation**
1. Go to the usearch download homepage of the 32-bit version (www.drive5.com/usearch/download.html).
2. Select the version "USEARCH v.10.0.240", select "Linux", register your email adress.
3. USEARCH will be sent to you by mail as an executable file.
3. Copy the usearch file into the order where you want to perform your analysis.
4. Rename the usearch file to "US_10_240" by typing "mv usearch10.0.240_i86osx32 US_10_240".
5. Run the command ""chmod +x US_10_240" and type in your password to make the file executable

**Vsearch installation**
For details, go to the vsearch homepage (https://github.com/torognes/vsearch).

```wget https://github.com/torognes/vsearch/archive/v2.8.0.tar.gz
tar xzf v2.8.0.tar.gz
cd vsearch-2.8.0
./autogen.sh
./configure
make
make install  # as root or sudo make install
wget https://github.com/torognes/vsearch/archive/v2.8.0.tar.gz
tar xzf v2.8.0.tar.gz
cd vsearch-2.8.0
./autogen.sh
./configure
make
make install  # as root or sudo make install
```
# FASPA script collection
First, go into the folder where you want to perform your analysis. Then, download the FASPA1 scripts to Linux by writing
```
wget https://github.com/StefanPfeiffer80/FASPA.github.io/
tar xzf v2.8.0.tar.gz
```

# Amplicon Processing using the bash scripts *FASP_preprocess.sh, FASP_unoise.sh, FASP_uparse.sh*
When all needed programs are installed and they are also at their place, amplicon processing with FASPA is pretty easy.
1. Go to your folder where you want to perform your analysis.
2. List your files by typing "ls" or "ll". Your folder should contain the following :
```
FASP_preprocess.sh
FASP_preprocess_vsearch.sh
FASP_unoise.sh
the usearch executable "US_10_240"
```
Put the bash scripts in the folder where your fastq files are. Your folder should then look somehow like this (+ of course most likely a higher number of .fastq files).
<p>
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/FASP_folder.png" width="620" height="100" />
</p>

# Amplicon Processing using the bash scripts *FASP_preprocess.sh, FASP_unoise.sh, FASP_uparse.sh -setting parameters*
Now that all files are in place, we have to configure the files that have to be configured. In this tutorial, we assume that you have a huge number of data, several hundred .fastq files. Also, we assume that we want to denoise raw reads("FASP_unoise.sh") instead of clustering OTUs ("FASP_uparse.sh"). For this reason this tutorial will take 

**FASP_Preprocessing.sh** is a bash script for the preprocessing of raw fastq files based on the programs USEARCH v10.240 and VSEARCH v2.80. Run the script by typing:

```
bash FASP_preprocess_vsearch.sh
```
After a while, the script will ask you if you have defined the length of the primers and the expected sequence length.
Thus we open the **"FASP_precrocess_vsearch.sh"** script using a text editor (just right click and choose "open with text editor" / if you are using the terminal only you can open by typing:"nano FASP_preprocessing.sh").
Now you have to look at the positions that are marked whether with "XX" substitute them with the length of your primers. "XXX" have to be replaced with the minimum and maximum expected length of your sequences. See the screenshot for a better understanding:

**IMPORTANT!!!!!!! Positions marked with XX need to be adjusted according to the users need!!!!!!!**
<p>
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/preprocess_selection.png" width="1000" height="80" />
</p>

After this is done, save your script and simply run the script as described above.

What happens by running **FASP_preprocess.sh** or **FASP_preprocess_vsearch.sh**?

a. Forward and reverse paired-end reads are merged and then all merged reads are put into one single .fastq file. The output file is named **raw.fq** -> more information http://drive5.com/usearch/manual/merge_pair.html. Further, an info file on the single merged file is created, named **raw_info.txt**

b. Expected Error rates (EE values) are created which indicates the probability if a particular base is right or wrong. The output file is named **qualrawfq.txt**. If you open the script in a text editor, you can change the calue to whether 0.5 (more stringent) or 2.0 (less stringent). For more information  on http://drive5.com/usearch/manual/exp_errs.html 

c. Trimming of primers, overhangs and quality filtering of the reads based on estimated error rates. The output file if we execute **FASP_preprocess_vsearch** is **vsearchfilteredstripped.fa**, while for **FASP_preprocess.sh** there are several output files: **strippedraw.fq** following trimming, and **filteredstripped.fa** following the subsequent quality filtering.  

d. Generation of a fasta file containing only unique sequences, the output file is **uniques.fa**

**FASP_unoise.sh** or **FASP_uparse.sh**
In FASPA, you can choose whether you want to denoise your filtered and trimmed raw reads or if you want to cluster OTUs at97% similarity level. 
    Unoise3                                            |       Uparse
-------------------------------------------------------| --------------------------------
Distinguish correct biological sequence from noisy read| Cluster similar reads into OTUs
no treshold, denoised raw reads                        | 97% similarity treshold

a. Using the uniques.fa file as an input, the unoise3 algorithm will create denoised raw reads, so called Zero-radius OTUs (ZOTUs). Further, a raw OTU table is created that is quality checked and 

b. The raw OTU table is further processed using the USEARCH's UNCROSS algorithm to get rid off wrongly assigned OTUs through cross-talk  More *information see: http://drive5.com/usearch/manual/crosstalk.html*
c. Taxonomic assignment of the OTUs is done using the SINTAX algorithm (paper) implemented in USEARCH and the rdp_16s_v16.fa database. The SINTAX algorithm uses k-mer similarity (https://en.wikipedia.org/wiki/K-mer) to identify the highest taxonomic ranks and provides an output table with bootstrap confidence values for all predicted taxonomic ranks.

 

Further, a phylogenetic tree in Newick-format (https://en.wikipedia.org/wiki/Newick_format is constructed via creation of a distance matrix and agglomerative clustering.
Output files of FASP_unoise.sh: -> zotus.fa

3.	FASP_uparse.sh
Using the uniques.fa file as an input, the uparse algorithm will create operational taxonomic units (OTUs) at a 97% sequence similarity treshold. 97% 16S rRNA sequence similarity was for a long time treated as a marker for the establishment of species. However, it should be considered that this will not be the case for many bacterial species, including most enterobacteriaceae., . Using OTUs defined by clustering from unique sequences, FASPA creates a OTU table, likewise corrected by applying the UNCROSS algorithm, as it was done in FASP_unoise.sh.

Output files of FASP_uparse.sh: -> zotus.fa
or OTU clustering (using UPARSE), both part of Usearch v.10.240

# Statistical analysis using R-studio and Rhea
Files which were generated by the workflow can be further analyzed using R, the most commonly used statistical language. While it is possible 
For a first statistical overview, I recommend the R-script collection Rhea (Lagkouvardos et al. 2016). Rhea can be downloaded from this link (LINK). Here, a quick tutorial is provided how you can implement your FASP generated files into Rhea.


Transparency: FASPA is a completely transparent workflow, advantegous for the user, also gives respect to the used programs.
4.	SINTAX algorithm for the taxonomic assignment and classification of OTUs or ASVs.
5.	Custom perl and R scripts that allow downstream analysis of the generated dataset i downstream analysis to utilize the statistical features of the popular tools QIIMEX the  R-based package phyloseqX and the Rhea-pipeline (L) which allows 
