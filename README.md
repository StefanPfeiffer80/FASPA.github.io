<p>
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/FASP_logo.png" width="800" height="150" />
</p>

        
**<p align="center">
    What is FASPA?
    </p>**
    
- FASPA is a workflow for analysing Illumina paired-end sequence data.    
- FASPA is a collection of shell bash scripts, perl scripts and R scripts and relies on state of the art programs used in sequence  
- processing, USEARCH and VSEARCH. FASPA output files are directly in downstream analysis of the data using the phyloseq package in R, the Rhea script collection, or the QIIME software package.   
- FASPA is distributed without warranty.  
- For using FASPA scripts, please cite: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses. DOI:kommt noch 
- If you have any questions, critics or suggestions, please write me at microbiawesome@gmail.com.

**<p align="center">
    FASPA - Fast Amplicon Sequence Processing and Analysis for MiSeq paired end sequence data
        </p>**
# Introduction
High-throughput sequencing of the 16S rRNA marker gene is the current benchmark in the characterization of bacterial microbial communities from virtually all environments.  
Here, I present **FASPA - Fast Amplicon Sequence Processing and Analysis**, an amplicon processing workflow that is easy to use and allows in depth-analysis of microbial communities. 
FASPA was created to address the need for a pipeline that can be used without prior bioinformatic knowledge, provides the user full transparency on which programs were applied and which parameters were used and at the same time allows to adjust and expand the pipeline according to the users needs.  
FASPA is based on executable bash scripts, perl scripts and allows very fast processing (less than 30 minutes for an avergae dataset of ~100 samples on an average speed laptop with 4 GB RAM), and is thus very well suited for the high-throughput analysis of 16S amplicon MiSeq data.
FASPA utilizes state of the art versions of much used bioinformatics tools USEARCH (Edgar 2010) and VSEARCH (Rognes et al. 2016) for amplicon processing and gives the user full transparency on how the tools are applied. 
FASPA also supports the integration of the processed amplicon data in various popular analysis tools, such as QIIME (Caporaso et al. 2010) or the R-based (R-core team, 2008) package collection phyloseq (McMurdie and Holmes, 2013). To gain a fast and detailed overview of the data, FASPA includes scripts that format the processed sequence data for direct application into the Rhea pipeline (Lagkouvardos et al. 2016).
In a nutshell, FASPA manages a precarious balance by being very fast, state of the art, having low CPU requirements, being transparent and at the same time easy to use.   

**Software needed**  
The FASPA script collection can be downloaded here. As FASPA relies on third party software workflow, it is necessary that you also install USEARCH v.10.240 (Edgar 2010) and optionally VSEARCH v.2.80 (or higher). 
USEARCH by Robert Edgar is a program that harbors a collection of functions and algorithms to efficiently, fast and accurately transform raw amplicon reads into an OTU table for downstream analysis. Many popular pipelines such as QIIME (up to v.1.9.), LotuS, IMNGS among others rely on different versions of USEARCH and other programs. USEARCH however, is free of charge only in the 32 bit version, which holds a 4GB memory cap. While this is not a problem for most datasets (depending on the sample type between 50 and 100 samples can be processed with the 32bit version), larger datasets will not be processed. USEARCH can be downloaded as a single executable file (www.drive5.com/usearch/download.html), see the tutorial for detailed instructions. To analyze larger datasets, FASPA includes VSEARCH by Torbjørn Rognes (Rognes et al. 2016), which was designed as an open source alternative to USEARCH (both USEARCH and VSEARCH are written in C++). For the FASPA workflow, it is evident that you don't use any VSEARCH version prior to v2.8.0.

# Which publications have to be cited when using FASPA?

For using FASPA scripts, please cite: Pfeiffer, S. (2018) FASPA - Fast Amplicon Sequence Processing and Analyses. (DOI:kommt noch)  

FASPA is transparent on which third party software it relies on.  I think it is important for the user to know which and how programs in the pipeline are applied, and it is fair for the developer of the third party software to be credited accordingly. 

Thus, in all cases when you use FASPA, cite also the USEARCH software, and also its implemented algorithms UNCROSS and SINTAX.
- Edgar,RC (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461.
doi: 10.1093/bioinformatics/btq461
- Edgar, R.C. (2016), UNCROSS: Filtering of high-frequency cross-talk in 16S amplicon reads. doi: http://dx.doi.org/10.1101/088666
- Edgar, R.C. (2016), SINTAX, a simple non-Bayesian taxonomy classifier for 16S and ITS sequences, http://dx.doi.org/10.1101/074161.

If you choose to preprocess your raw reads using *FASP_preprocessing_v1.sh* (the setup for large datasets using VSEARCH), cite:
- Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584

If you choose to denoise your sequencesby calling *FASP_unoise.sh*, cite:
- Edgar, R.C. (2016), UNOISE2: Improved error-correction for Illumina 16S and ITS amplicon reads.http://dx.doi.org/10.1101/081257

If you choose to cluster your sequences into OTUs by calling *FASP_uparse.sh*, cite:
- Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods https://www.nature.com/articles/nmeth.2604

If you choose to perform downstream analysis in **R** and **RStudio**, cite:
- R-core team R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org.
- RStudio Team (2015). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/.
The **Rhea** script collection:
- Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 https://doi.org/10.7717/peerj.2836
**Phyloseq**: 
McMurdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217
**RDPutils**
John Quensen (2018). RDPutils: R Utilities for Processing RDPTool Output.

# FASPA Tutorial
# 1. Before we start
This tutorial intends to show an user with verly low bioinformatic experience how to process and analyze 16S amplicon sequencing data on her/his personal computer. As FASPA is based not on Micosoft Windows, but on the Linux operating system, I recommend the installation of a virtual box on your computer, where you can then run a Linux distribution such as Ubuntu. 
Follow this link to set up your system: https://www.wikihow.com/Install-Ubuntu-on-VirtualBox

**Download and extract the FASPA script collection**
First, go into the folder where you want to perform your analysis. This should be the same folder where your fastq-files are. Then, download the FASP_1.0.tar.gz file to Ubuntu by writing:
```
wget https://github.com/StefanPfeiffer80/FASPA.github.io/FASP_1.0.tar.gz
tar xzf FASP_1.0.tar.gz
```

**Change the primers.fa file according to your needs**
You will need a file in fasta format that contains the names and sequences of your primer pair. Open the primers.fa file that you downloaded through FASP_1.0.tar.gz and exchange the primer names and the primer sequences. The forward primer should be always the upper two lines. Then save the file without changing the name.
<p align="center">
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/prmers.fa.png" width="300" height="160" />
</p>
**USEARCH installation**
1. Go to the USEARCH homepage where you can download the 32-bit version (www.drive5.com/usearch/download.html).
2. Select the version "USEARCH v.10.0.240", select "Linux", register your email adress.
3. USEARCH will be sent to you by mail as an executable file.
3. Copy the USEARCH file into the order where you want to perform your analysis.
4. Rename the USEARCH file to "US_10_240" by typing "mv usearch10.0.240_i86osx32 US_10_240" or just rename the file in the GUI.
5. Run the command ""chmod +x US_10_240" and type in your password to make the file executable

**VSEARCH installation**
1. Open a terminal and copy paste the text in the box. Keep in mind that you need admin rights to install VSEARCH.  VSEARCH is updated quite frequently. Thus I recommend to follow the installation instructions at the VSEARCH homepage (https://github.com/torognes/vsearch).

**16S database for taxonomic assignment**
FASPA by default uses the RDP_16S_v16 which is also a recommendation for the SINTAX classifier used in FASPA. To download the training set, you can go to the USEARCH homepage: https://www.drive5.com/usearch/manual/sintax_downloads.html.  
```
wget https://www.drive5.com/sintax/rdp_16s_v16.fa.gz
gunzip rdp_16s_v16.fa.gz
```

# 2. FASPA workflow tutorial
**Amplicon Processing using the bash scripts *FASP_preprocess_v1.sh, FASP_preprocess_v2.sh, FASP_unoise.sh, FASP_uparse.sh***  

When all needed programs are installed and they are also at their place, amplicon processing with FASPA can finally start.
1. Go to your folder where you want to perform your analysis.
2. List your files by typing "ls" or "ll". Your folder should then look somehow like this (there will be different and of course I assume a higher number of .fastq files). Bild muss mit filenames updated werden!!!
<p align="center">
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/FASP_folder.png" width="620" height="100" />
</p>

Now that all files are in place, we have to configure the files that have to be configured. 

**FASP_preprocess_v1.sh** and **FASP_preprocess_v2.sh**
The preprocessing step prepares the raw amplicon reads for subsequent OTU clustering or denoising.
Preprocessing includes:
- merging of paired reads
- quality filtering
- primer removal
- length trimming
- preparing a set of unique sequences

For small and medium sized 16S amplicon libraries of up to 100 samples, the script *FASP_preprocess_v1.sh* can be called. The advantage is that vsearch does not need to be installed, as *FASP_preprocess_v1.sh* uses the 32 bit version of USEARCH which has a 4GB memory cap. Thus, *FASP_preprocess_v1.sh* should be only used if you have a relatively small dataset.  
For running *FASP_preprocess_v1.sh*, you have to set to tell FASPA the length of the primers you are using primer. In the example below we use the primer pair by Nossa et al. (2010), wheras the forward primer (-l) consists of 20 nucleotides, and the reverse primer (-r)  of 19. 
```
bash FASP_preprocess_v1.sh -l 19 -r 20
```

If we assume that you have a huge number of data, several hundred .fastq files, we use *FASP_preprocess_v2.sh*, which beside USEARCH, applies also the open source software VSEARCH for the sequence processing. In order to run bash *FASP_preprocess_v2.sh*, we have to set several parameters.  
These include
- -l the length of the forward primer
- -r the length of the reverse primer
- -m the maximum expected sequence length
- -s the minimum expected sequence length

In the example below we use again the primer pair by Nossa et al. (2010), which confines an expected amplicon length of ~450 bp. Both primers have a length of 19 nucleotides, and thus the the expected length of the sequence is around 410 bp
```
bash FASP_preprocess_v2.sh -l 19 -r 19 -m 450 -s 390
```
Now the script starts running. We see that forward and reserve reads are merged, and we see the percentage of reads that are merging. All files will be stored together in one file, named **raw.fq**. Next, FASPA extracts randomly a subset of your reads (by default 100). This is done to check at which position of your sequences the actual primer sequences are found / or not found. At this point the script stops and asks you if you checked the position and size of your primers and the expected length of your amplicons:
<p align="center">
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/check.png" width="620" height="100" />
</p>
To answer these questions, we look in our folder and see that a new file appeared, named: *primer_positions.txt*. Open the file with a text editor (in ubuntu by right click, in the terminal by typing "nano primer_positions.txt"). Here you see now a table with four columns:

<p align="center">
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/primrcheck.png" width="620" height="100" />
</p>

In the left column there is the identifier of the merged readpair, e.g.Stefan5.43832. In the next column, you find the starting position of the primer and in the third column the end position. According to our primer length that we know from the primers.fa file, there should be in our case 19 bases difference between the starting positions and the end position. In the fourth column we see that the primer which was found at positions 1-19 aligned to the + strand / R1 read and the primer which was found at positions >400 aligned to the - strand / R2 read. We can see that not all reads have exactly the same length. That is because we look at a bacterial community, with reads from different microorganisms. Thus, differences in read length are explained via deletions or insertions that characterize different bacterial lineages. For this reason, when you start the script, you should choose minimum and maximum values in an at least +/- 50bp range around the expected length.  
Another case that will likely occur when you look into the *primer_positions.txt* list is that for one read only one primer was found. In this case, the read is too short that both primers could be detected.

After we checked that the primers are actually there and found in the expected length, we can choose 1) and press the enter key. However, if you realized that the provided primer length was wrong or that the expected size of amplicons is lower or higher, you should stop the script here by choosing option 2, and restart the script with the corrected parameters.

The script continues. First, primers are stripped from the sequences and sequences are filtered for the selceted length. Next, reads that don't reach the minimum quality requirements are filtered out.". Last, unique sequences are extracted and saved as *uniques.fa*. This is done to significantly reduce the datasize for the following denoising or clustering algortihms (Abundance data on how often each read was found is still integrated in the output file). When the script is finished successfully, FASPA tells you that you can continue the amplicon processing by whether denoising of the cleaned raw reads or clustering of OTUs.

*Output files*
- *raw.fq*          All paired reads merged and stored in one file. More information: http://drive5.com/usearch/manual/merge_pair.html
- *raw_info.txt*    Provides information on the length distribution of reads (min, max, median, average) and the mean estimated error rates
- *qualrawfq.txt*   EE values (=sum of error probabilities) are created to indicate the probability contains errors. The default set in the script is set to 1.0, which corresponds to zero errors. You can manually change the parameter when you open the script in a text editor. Lower values e.g. 0.5 are more stringent. For more information: http://drive5.com/usearch/manual/exp_errs.html. 
- *strippedraw.fq*  The *raw.fq* file with the primer sequences removed.
- *filtered.fa*     The *strippedraw.fq* file with low quality reads removed. Note that this is a fasta file, not a fastq file.
- *uniques.fa*      This file contains the uniques sequences of *filtered.fa*.

**FASP_unoise.sh** or **FASP_uparse.sh**
In FASPA, you can choose whether you want to denoise your preprocessed raw reads or if you want to cluster OTUs at 97% sequence similarity level.
USEARCH includes the UPARSE algorithm (Edgar 2013) to cluster OTUs, which showed improved accuracy in OTU assignment towards other commonly used clustering algorithms and was already cited several thousand times. For details on the clustering algorithm, see here: https://www.drive5.com/usearch/manual/uparseotu_algo.html. In version 9, USEARCH implemented UNOISE (Edgar and Flyvbjerg 2015, Edgar 2016), an algorithm for denoising of raw sequences, which actually means that the genetic variation of sequences is analyzed to find out what causes the sequence variation; whether real sequence differences or sequencing errors. For details see https://www.drive5.com/usearch/manual/unoise_algo.html.  
Today, denoising of raw amplicon reads becomes more popular especially in hindsight of the known biases that go together with OTU clustering using a 97% sequencing similarity cutoff to differentiate between species. For more information, look at the reviews !!(XXXX)). Under perfect circumstances denoising will always yield more OTUs than clustering at the 97% species level, as different strains of one species should be identified by denoising (which is not the case in reality often). To figure this out, I recommend to have one mock community of known bacterial strains in your samples to analyze, and compare the results of denoising with OTU clustering.
In this tutorial, I expl

**FASP_unoise.sh** and **FASP_uparse.sh**
In order to run the bash-script *FASP_unoise.sh* or *FASP_uparse.sh*, we have to set one parameter:
- -i the minimum length of your ZOTUs/OTUs

```
bash FASP_unoise.sh -i 350
#OR
bash FASP_uparse.sh -i 350
```

Using the *uniques.fa* file as an input, the UNOISE3 algorithm will create denoised raw reads, so called Zero-radius OTUs (ZOTUs). 
Then the script stops again and asks you if you have checked the expected length of your ZOTUs or OTUs? 
This is most important when you were running *FASP_preprocess_v1.sh*, as there was not length trimming done beforehand and your file might contain relatively short ZOTUs/OTUs. If everything is alrgight, choose option 1) and press the enter key. 

Next FASPA creates an OTU table (or ZOTU table, but for simplicity we will use OTU for both ZOTU and OTU from now on) out of the OTUs and the raw reads. The OTU table is a ext file which can be opended in any spreadsheet program and also coincides to the format of QIIME classic OTU tables. Sample IDs are displayed in columns and OTUs are displayed in rows.
Next, the raw OTU table is further processed using the USEARCH's UNCROSS (Edgar 2016) algorithm to get rid off wrongly assigned OTUs through cross-talk. For more information crosstalk see: http://drive5.com/usearch/manual/crosstalk.html. Taxonomic assignment of the OTUs is done using the SINTAX algorithm (Edgar 2016) and the rdp_16s_v16.fa database. The SINTAX algorithm uses k-mer similarity (https://en.wikipedia.org/wiki/K-mer) to identify the highest taxonomic ranks and provides an output table with bootstrap confidence values for all predicted taxonomic ranks. If you want to use another database, adjust the input file for the taxonomic databae "-db" at #7 of the FASP_unoise.sh/FASP_uparse.sh script. Also you can adjust the cutoff value "-sintax_cutoff" which is set to 0.5 by default (corresponds to 50% bootstrap support).

Maybe a picture here?

In FASPA, SINTAX assigned OTUs are filtered from the taxonomy file for OTUs that are not occurring in the OTU table (e.g. OTUs that were filtered out following the UNCROSS algorithm) to avoid cover inequality of the OTU table and the taxonomic data. Next, a phylogenetic tree in Newick format (https://en.wikipedia.org/wiki/Newick_format is constructed via creation of a distance matrix and agglomerative clustering, that will be later used to calculate UNIFRAC distances for diversity analyses.
Last, two scripts are called that create OTU tables with added taxonomic information. One table, *otutab_otus_greengenes.txt* contains taxonomic information in the greengenes syntax, the other one, *otutab_otus_SILVA.txt*, contains taxonomic information using the SILVA syntax. These files might be useful if you want to perform downstream analysis in QIIME or if you want to predict hypothetical metagenomic information with PICRUSt (Langille et al. 2013).

*Output files*      *UN* as in *utUN.txt* relates to *FASP_unoise.sh*, *UP* as on *outUP.txt* relates to*FASP_uparse.sh*
- *ZOTUS.fa* or *OTUS.fa*         Fasta sequences of denoised reads (Zotus.fa) or clustered OTUs (Otus.fa)
- *otus_sorted_UN.fasta*    OTUs filtered by minimum length set by parameter -i
- *otutab_UN_raw.txt*       Raw OTU table
- *otutab_UN_uncrossed.txt* OTU table corrected for sequencing errors by crosstalk-filtration.
- *otutab_UN_uncrossed_stats.txt*               Statistics of the OTU table *otutab_UN_uncrossed.txt*. Contains information such as total number of reads, number of OTUs and core OTUs (OTUs that occur in 100%, 90% and 50% of samples)
- *sintaxzotusrdp_UN.txt* Taxonomic assignment of OTUs (*otus_soretd_UN_fasta*) using the rdp_16s_v16.fa database (default), including bootstrap values generated.
- *SINTAX_OTUS_RDP_UN_FILT.txt*    Adjusted output of *sintazotusrdp.txt* to match the number of OTUs in the OTU table, after performing the *FASP_tax_filtered.pl* command.
- *Tree_UN.tree*    A tree file in the Newick format that can be used in downstream analyis.
- *otutab_otus_greengenes.txt*   The otu table *otutab_UN_uncrossed.txt* with taxonomic information added using the greengenes syntax.
- *otutab_otus_SILVA.txt*   The otu table *otutab_UN_uncrossed.txt* with taxonomic information added using the SILVA syntax.

# 3. Removal of contaminant OTUs
The addition of negative controls for PCR (and also of negative extraction controls to determine contaminants from extraction kits or sample processing) is crucially important in amplicon library preparation. When it comes to the analyis of the sequencing data, we want to filter out OTUs that derive from the negative controls and contaminate our real biological samples. Although the removal of these contaminant OTUs could be easily automated and included in the FASP scripts, I recommend to investigate your OTU table via a spreedsheet program such as Microsoft Excel or LibreOffice Calc. The reason for this is that by automated filtration you might lose OTUs that may appear 10000 times in particular samples but just once or twice in your control. Knowing that especially highly abundant sample OTUs can appear as false positives in the controls via e.g. crosstalk, it is definitely a good choice to investigate your dataset manually. In FASPA, although these typically wrongly assigned crosstalk-OTUs are filtered out mostly beforehand thanks to the UNCROSS algorithm (check this out by comparing *otutab_UN_raw.txt* with *otutab_UN_uncrossed.txt*), only a manual check-up will give you certainity. When you identified your contaminant OTUs, simply delete them from all samples and then also the now empty control samples. Save the table with a new name, e.g. *otutab_UN_uncrossed_filt.txt*.  
Then run the script *FASP_tax_filtered.pl*

```
perl ./FASP_tax_filtered.pl otutab_UN_uncrossed_filt.txt SINTAX_OTUS_RDP_FILT.txt SINTAX_OTUS_RDP_FILT2.txt
```
By this script, OTUs that were removed by taking out contaminants from the negative control, are also taken out from the taxonomy file.
(output = *SINTAX_OTUS_RDP_FILT2*).

# 4. Statistical analysis using R-studio and Rhea
Files which were generated by the workflow can be further analyzed using the widely used statistical language R (R Development Core Team, 2008). I recommend to perform the downstream analysis in the user-friendly Rstudio (RStudio Team 2015). 
Phyloseq (McMurdie and Holmes, 2013) enjoys great popularity as it combines in depth statistical analyses of amplicon sequencing data with high-quality graphs. To gain a statistical overview however via a stanardized statistics pipeline, I suggest to use the R-script collection Rhea (Lagkouvardos et al. 2016). 
A common problem in bioinformatics is how to combine different bioinformatic analyis tools due to differences in file types, the data structure and the syntax used.  
*FASPA* contains several R-scripts that allow transformation of the output files from **FASP_uparse.sh** and **FASP_unoise.sh** into a single *phyloseq-object*. From there on, the script **Physeq_to_Rhea.R** further transforms the *phyloseq-object* into a *OTU table with taxonomy* which is already formatted according to the needs of the Rhea pipeline. Additionally, we need to install the package tidyr (select from your closest cran repository).  
In the following lines, I will provide an exact guideline on how out of your output files from **FASP_uparse.sh** amd **FASP_unoise.sh** a *phyloseq-object* and a *Rhea-input-file* will be created.  
# R-Tutorial 
**Software needed**
Go to the R-studio page: https://www.rstudio.com/products/rstudio/download/. Choose your operating system, download the newest version, then follow the installation instructions.
Then start R-studio and install phyloseq (see the detailed instructions here: https://joey711.github.io/phyloseq/install.html). 
Rhea can be downloaded from this link (https://lagkouvardos.github.io/Rhea/). With the downloaded files you will also find a .pdf file with detailed instructions on how to apply each of the Rhea scripts. We also need the package RDPutils (Quensen 2018), whcih can be installed following the instruction on this site (https://rdrr.io/github/jfq3/RDPutils/). Additionally, we need to install the package tidyr (Wickham and Henry 2018; select from your closest cran repository). Finally, download the FASPA *Physeq_to_Rhea.R* script file from this repository.
Create the file location (In Windows or Linux) where you want to perform your analysis. e.g. C:\RheaTest\.
In this folder you have to prepare six empty folders, named exactly as shown in the picture below.

<p align="center">
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/" width="620" height="100" />
</p>

**FASPA output files needed**
For simplicity, the tutorial deals with the output files of *FASP_unoise.sh*. 
- *Tree_UN.tree*
- *SINTAX_OTUS_RDP_UN_FILT.txt*
- *otutab_UN_uncrossed.txt*
- *mapping_file.tab*  The mapping file contains all metainformation of your samples, including categorical and numerical variables. You can find an example mapping file here.

Copy all output files in the 1.Normalization folder. Also, copy the files *Tree_UN.tree* and *mapping_file.tab* in the "Beta-Diversity" folder. Also deposit one copy of *mapping_file.tab* into the "5.Group-Comparisons" folder

**Preparation of FASPA output files for phyloseq and Rhea**
First we will transform the input files into a *phyloseq-object*, which we name *physeq* (see the code block below). As the phyloseq package is loaded you can perform various analyses already.
You have to adjust the names in the code block below if you used different names for your FASP output files (e.g. after filtering contaminants).
Run the code below:

```
setwd("C:/RheaTest/Rhea/1.Normalization")
TAX<-RDPutils::import_sintax_file("SINTAX_OTUS_RDP_UN_FILT.txt")
otumat<-read.csv("otutab_UN_uncrossed.txt",sep="\t", header = TRUE, row.names=1)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
physeq = phyloseq(OTU, TAX)
sampledata = import_qiime_sample_data("mapping_file.tab")
physeq = phyloseq(OTU,TAX,sampledata)
OTU_TAX_Rhea <- cbind.data.frame(otu_table(physeq),tax_table(physeq))
```

Then, we open the script *Physeq_to_Rhea.R* in RStudio.
In the shell, we look at the structure of the OTU+taxonomy table that was extracted out of the phyloseq object by typing "str(OTU_TAX_Rhea)" (see the code block below). The output will give you the number of observations (rows = OTUs) and variables (columns = samples + taxonomic ranks). 

Now, look at the *Physeq_to_Rhea* script in the editor: Here you have to adjust the positions of your taxonomic columns, which will be the last 6 variables/columns (e.g. if *str(OTU_TAX_Rhea)* showed 33 variables, this will be the variables "28:33". Accordingly, you have to adapt the number of samples (see the line that starts with R4!) In our case the positions will be 1:27.

```
Physeq_to_Rhea<-function(OTUTAX){
  R1<-OTUTAX[,28:33#<-Select columns which contain the (normally 6) taxonomic levels = Starts with Number of samples+2 
             ]%>%unite(taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";",remove=TRUE) #<-Adjust 
##code not shown##
        R4 <- gsub(";","/t", R3[,1:27])#<-Number of samples counting from 1 + 1
##code not shown##
```

After you adjusted *Physeq_to_Rhea* script, copy & paste the content in the shell and press the Enter key!
Then type:
```
Physe_to_Rhea(OTU_TAX_Rhea)
```
This will create an output file *OTU_table_w_taxonomy.txt* which can be directly applied in the Rhea pipeline, starting with the *Normalization* script (see the instructions in the Rhea script files).

# 6. Changing parameters
FASPA script files can be adjusted to your personal needs and preferences. 
<p>
    <img src="https://github.com/StefanPfeiffer80/FASPA.github.io/blob/master/pictures/preprocess_selection.png" width="1000" height="80" />
</p>


# 7. Citations
- Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, et al. (2010) QIIME allows analysis of high-throughput community sequencing data,  Nature Methods, doi:10.1038/nmeth.f.303  
- Edgar RC. (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461.
doi: 10.1093/bioinformatics/btq461.
- Edgar RC. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods, doi:10.1038/nmeth.2604.
- Edgar RC. (2016), UNCROSS: Filtering of high-frequency cross-talk in 16S amplicon reads. doi: doi:10.1101/088666.
- Edgar RC. (2016), SINTAX, a simple non-Bayesian taxonomy classifier for 16S and ITS sequences, doi:10.1101/074161.
- Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584, https://github.com/torognes/vsearch
- Edgar RC. (2016), UNOISE2: Improved error-correction for Illumina 16S and ITS amplicon reads.doi:10.1101/081257
- Lagkouvardos I, Fischer S, Kumar N, Clavel T. (2017) Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons. PeerJ 5:e2836 doi:10.7717/peerj.2836
- Murdie and Holmes (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 8(4):e61217, https://joey711.github.io/phyloseq/
- Quensen J. (2018). RDPutils: R Utilities for Processing RDPTool
  Output. R package version 1.4.1. https://rdrr.io/github/jfq3/RDPutils/
- R-core team R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org.
- RStudio Team (2015). RStudio: Integrated Development for R. RStudio, Inc., Boston, MA URL http://www.rstudio.com/
- Wickham and Henry (2018)tidyr: Easily Tidy Data with 'spread()' and 'gather()' Functions. R package version 0.8.0.
  https://CRAN.R-project.org/package=tidyr



 


3.	FASP_uparse.sh
Using the uniques.fa file as an input, the uparse algorithm will create operational taxonomic units (OTUs) at a 97% sequence similarity treshold. 97% 16S rRNA sequence similarity was for a long time treated as a marker for the establishment of species. However, it should be considered that this will not be the case for many bacterial species, including most enterobacteriaceae., . Using OTUs defined by clustering from unique sequences, FASPA creates a OTU table, likewise corrected by applying the UNCROSS algorithm, as it was done in FASP_unoise.sh.

Output files of FASP_uparse.sh: -> zotus.fa
or OTU clustering (using UPARSE), both part of Usearch v.10.240
