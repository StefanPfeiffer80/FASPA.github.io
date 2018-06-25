
# FASPA - Fast Amplicon Sequence Processing and Analysis for MiSeq paired end sequence data

If you have any questions, critics or suggestions, please write me at microbiawesome@gmail.com.
FASPA is a workflow for analysing Illumina paired-end sequence data. 
FASPA is a collection of shell bash scripts, perl scripts and R scripts and relies on state of the art programs used in sequence processing, USEARCH and VSEARCH. FASPA output files are directly in downstream analysis of the data using the phyloseq package in R, the Rhea script collection, or the QIIME software package.   
FASPA is distributed without warranty.  
For using FASPA scripts, please cite: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses. DOI:kommt noch  

# 1. Introduction
High-throughput sequencing of the 16S rRNA marker gene is the current benchmark in the characterization of bacterial microbial communities from virtually all environments.  
Here, I present FASPA -Fast Amplicon Sequence Processing and Analysis, an amplicon processing workflow that is easy to use and allows in depth-analysis of microbial communities. 
FASPA addresses the need for a transparent pipeline that based on executable bash scripts, perl scripts and  which allows very fast processing (less than 30 minutes on an average speed laptop with 4 GB RAM).
FASPA utilizes the well-known and state of the art bioinformatics tools USEARCH (Edgar 2010) and VSEARCH (Rognes et al. 2016) and gives full transparency on how the tools are applied. 
FASPA also supports the integration of the processed amplicon data in various popular analysis tools, such as QIIME (Caporaso et al. 2010) or the R-based (R-core team, 2008) package collection phyloseq (McMurdie and Holmes, 2013). To gain an overview of the amplicon data (recommended), FASPA includes scripts that format the processed sequence data for direct application of the Rhea pipeline (Lagkouvardo et al. 2016) in R (R-core team, 2008).
In a nutshell, FASPA manages precarious balance by being very fast, applies state of the art bioinformatic tools, having low CPU requirements, offers full transparency for the user by being at the same time easy to use. 

**Software needed**
For the first part of the FASPA workflow, the amplicon processing, FASPA uses the programs USEARCH v.10.240 and optionally VSEARCH v.2.80. 
USEARCH by Robert Edgar (Edgar 2010) is a collection of functions and algorithms to efficiently, fast and accurately transform raw amplicon reads into an OTU table for downstream analysis. Usearch can be downloaded as a single executable file (www.drive5.com/usearch/download.html). USEARCH includes the UPARSE algorithm (Edgar 2013) to cluster OTUs, which showed improved accuracy in OTU assignment towards other commonly used clustering algorithms and was already cited several thousand times. For detailson the clustering algorithm, see here: https://www.drive5.com/usearch/manual/uparseotu_algo.html. In version 9, USEARCH implemented UNOISE (Edgar and Flyvbjerg 2015, Edgar 2016), an algorithm for denoising of raw sequences, which actually means that the genetic variation of sequences is analyzed to find out what causes the sequence variation; whther real sequence differences or sequencing errors. For details see https://www.drive5.com/usearch/manual/unoise_algo.html.  
Today, denoising of raw amplicon reads becomes more popular especially in hindsight of the known biases that go together with OTU clustering using a 97% sequencing similarity cutoff to differentiate between species. For more information, look at the reviews XXXX  .
USEARCH however, is free of charge only in the 32 bit version, which holds a 4GB memory cap. While this is not a problem for most datasets (depending on the sample type between 50 and 100 samples can be processed with the 32bit version), larger datasets will not be processed.  For this reason, FASPA includes VSEARCH by Torbjørn Rognes (Rognes et al. 2016), which was designed as an open source alternative to USEARCH (both USEARCH and VSEARCH are written in C++). For the FASPA workflow, it is evident that you don't use any VSEARCH version prior to v2.8.0.

# 2. Which publications have to be cited

For using FASPA scripts, please cite: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses. DOI:kommt noch  

FASPA calls several functions of the USEARCH program, that have to be cited seperately.
In all cases, when you use FASPA, cite:
- Edgar,RC (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461.
doi: 10.1093/bioinformatics/btq461
- Edgar, R.C. (2016), UNCROSS: Filtering of high-frequency cross-talk in 16S amplicon reads. doi: http://dx.doi.org/10.1101/088666
UNOISE algorithm
- Edgar, R.C. (2016), SINTAX, a simple non-Bayesian taxonomy classifier for 16S and ITS sequences, http://dx.doi.org/10.1101/074161.

If you choose to preprocess your raw reads using FASP_preprocessing_v1.sh, cite:
Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

If you choose to denoise your sequencesby calling FASP_unoise.sh, cite:
Edgar, R.C. (2016), UNOISE2: Improved error-correction for Illumina 16S and ITS amplicon reads.http://dx.doi.org/10.1101/081257

If you choose to cluster your sequences into OTUs by calling FASP_uparse.sh, cite:
Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772,  dx.doi.org/10.1038/nmeth.2604].

If you choose to perform downstream analysis in R, cite:
R-studio citation -> RStudio Team (2015). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/.
The Rhea script collection ->
Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836
Phyloseq -> McMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217

# 3. Before we start
**Download and extract the FASPA script collection**
First, go into the folder where you want to perform your analysis. Then, download the FASPA1 scripts to Linux by writing
```
wget https://github.com/StefanPfeiffer80/FASPA.github.io/
tar xzf v2.8.0.tar.gz
```

**USEARCH installation**
1. Go to the usearch download homepage of the 32-bit version (www.drive5.com/usearch/download.html).
2. Select the version "USEARCH v.10.0.240", select "Linux", register your email adress.
3. USEARCH will be sent to you by mail as an executable file.
3. Copy the usearch file into the order where you want to perform your analysis.
4. Rename the usearch file to "US_10_240" by typing "mv usearch10.0.240_i86osx32 US_10_240".
5. Run the command ""chmod +x US_10_240" and type in your password to make the file executable

**VSEARCH installation**
1. Open a terminal and copy paste the text in the box. Keep in mind that you need admin rights to install VSEARCH.  For more information, go to the VSEARCH homepage (https://github.com/torognes/vsearch).

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
**16S database for taxonomic assignment**
FASPA by default uses the RDP_16S_v16 which is also a recommendation for the SINTAX classifier used in FASPA. To download the training set, you can go to the USEARCH homepage: https://www.drive5.com/usearch/manual/sintax_downloads.html.  
```
wget https://www.drive5.com/usearch/manual/sintax_downloads/rdp_16s_v16.fa.gz
tar xzf rdp_16s_v16.fa.gz
```

# 4. FASPA workflow tutorial
# Amplicon Processing using the bash scripts *FASP_preprocess_v1.sh, FASP_preprocess_v2.sh, FASP_unoise.sh, FASP_uparse.sh*
When all needed programs are installed and they are also at their place, amplicon processing with FASPA is pretty easy.
1. Go to your folder where you want to perform your analysis.
2. List your files by typing "ls" or "ll". Your folder should contain the following :
Put the bash scripts in the folder where your fastq files are. Your folder should then look somehow like this (+ of course most likely a higher number of .fastq files).
<p>
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/FASP_folder.png" width="620" height="100" />
</p>

# Amplicon Processing using the bash scripts *FASP_preprocess.sh, FASP_unoise.sh, FASP_uparse.sh -setting parameters*
Now that all files are in place, we have to configure the files that have to be configured. In this tutorial, we assume that you have a huge number of data, several hundred .fastq files. For this reason, we use FASP_preprocess_v1.sh, which beside USEARCH, applies also VSEARCH for the sequence processing. Also, we assume that we want to denoise raw reads("FASP_unoise.sh") instead of clustering OTUs ("FASP_uparse.sh"). For this reason this tutorial will take 

**FASP_Preprocessing.sh** is a bash script for the preprocessing of raw fastq files based on the programs USEARCH v10.240 and VSEARCH v2.80. Run the script by typing:

```
bash FASP_preprocess_2.sh
```
After a while, the script will ask you if you have defined the length of the primers and the expected sequence length.
Thus we open the **"FASP_precrocess_v1.sh"** script using a text editor (just right click and choose "open with text editor" / if you are using the terminal only you can open by typing:"nano FASP_preprocessin_v1.sh").
Now you have to look at the positions that are marked whether with "XX" substitute them with the length of your primers. "XXX" have to be replaced with the minimum and maximum expected length of your sequences. See the screenshot for a better understanding:

**IMPORTANT!!!!!!! Positions marked with XX need to be adjusted according to the users need!!!!!!!**
<p>
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/preprocess_selection.png" width="1000" height="80" />
</p>

After this is done, save your script and simply run the script as described above.

What happens by running **FASP_preprocess_v1.sh** or **FASP_preprocess_v2.sh**?

a. Forward and reverse paired-end reads are merged and then all merged reads are put into one single .fastq file. The output file is named **raw.fq** -> more information http://drive5.com/usearch/manual/merge_pair.html. Further, an info file on the single merged file is created, named **raw_info.txt**

b. Expected Error rates (EE values) are created which indicates the probability if a particular base is right or wrong. The output file is named **qualrawfq.txt**. If you open the script in a text editor, you can change the treshold on EE values whether to 0.5 (more stringent) or 2.0 (less stringent). For more information  on http://drive5.com/usearch/manual/exp_errs.html 

c. Trimming of primers, overhangs and quality filtering of the reads based on estimated error rates. The output file if we execute **FASP_preprocess_v1** is **vsearchfilteredstripped.fa**, while for **FASP_preprocess.sh** there are several output files: **strippedraw.fq** following trimming, and **filteredstripped.fa** following the subsequent quality filtering.  

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

# 4. Statistical analysis using R-studio and Rhea
Files which were generated by the workflow can be further analyzed using the widely used statistical language R (R Development Core Team, 2008), and I recommend performing the downstream analysis in the user-friendly Rstudio (RStudio Team 2015). 
Phyloseq (McMurdie and Holmes, 2013) enjoys great popularity as it allows as its output combines in depth statistical analyses of amplicon sequencing data with high-quality graphs.  
To gain a first statistical overview however, I recommend the R-script collection Rhea (Lagkouvardos et al. 2016). Rhea can be downloaded from this link (https://lagkouvardos.github.io/Rhea/). With the downloaded files you will also find a .pdf file with detailed instructions on how to use the files.
*FASPA* contains several R-scripts that allow transformation of the output files from **FASP_uparse.sh** and **FASP_unoise.sh** into a single *phyloseq-object*. From there on, the script **Physeq_to_Rhea.R** further transforms the *phyloseq-object* into a *OTU table with taxonomy* which is formatted according to the needs of the Rhea pipeline.  
In the following lines, I will provide an exact guideline how out of your FASPA output files a *phyloseq-object* and a *Rhea-input-file* will be created.  
# R-Tutorial
**Install R-studio**
Go to the R-studio page: https://www.rstudio.com/products/rstudio/download/. Choose your operating system, download the newest version, then follow the installation instructions.
Then start R-studio and install phyloseq (see the detailed instructions here: https://joey711.github.io/phyloseq/install.html). 
Additionally, we need to install the package tidyr. Select install packages 

**Preparing the files for downstream analysis**


Transparency: *FASPA* is a completely transparent workflow, advantegous for the user, also gives respect to the used programs.
4.	SINTAX algorithm for the taxonomic assignment and classification of OTUs or ASVs.
5.	Custom perl and R scripts that allow downstream analysis of the generated dataset i downstream analysis to utilize the statistical features of the popular tools QIIMEX the  R-based package phyloseqX and the Rhea-pipeline (L) which allows 
# 5. Citations
- Edgar,RC (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461.
doi: 10.1093/bioinformatics/btq461
- Edgar, R.C. (2016), UNCROSS: Filtering of high-frequency cross-talk in 16S amplicon reads. doi: http://dx.doi.org/10.1101/088666
UNOISE algorithm
- Edgar, R.C. (2016), SINTAX, a simple non-Bayesian taxonomy classifier for 16S and ITS sequences, http://dx.doi.org/10.1101/074161.

If you choose to preprocess your raw reads using FASP_preprocessing_v1.sh, cite:
Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

If you choose to denoise your sequencesby calling FASP_unoise.sh, cite:
Edgar, R.C. (2016), UNOISE2: Improved error-correction for Illumina 16S and ITS amplicon reads.http://dx.doi.org/10.1101/081257

If you choose to cluster your sequences into OTUs by calling FASP_uparse.sh, cite:
Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772,  dx.doi.org/10.1038/nmeth.2604].

R-core team R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org.

RStudio Team (2015). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/

Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836
 cMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217
QIIME allows analysis of high-throughput community sequencing data 
J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight; Nature Methods, 2010; doi:10.1038/nmeth.f.303 
