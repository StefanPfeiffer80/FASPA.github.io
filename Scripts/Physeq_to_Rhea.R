# Physeq to Rhea
# (C) Stefan Pfeiffer, 28/2/2018, last modified on 01/07/2018, all rights reserved.

setwd("D:/Users/R/Rhea") #<- set working directory where you want to write the Rhea output files

#load necessary packages
library(phyloseq)
library(tidyr)

#load your data
data("physeq1") #<-! change to the name of your phyloseq object!

#!!Adapt the function to your dataset, look at the arrows!!
# Function Physeq_to_Rhea 
# Transforms OTU table with taxonomy to Rhea format
Physeq_to_Rhea<-function(OTUTAX){
  R1<-OTUTAX[,17:22#<-Select columns which contain the taxonomic levels = last 6 columns
                    #Note: Rhea uses only 6 taxonomic ranks (down to "Genus").
                    #Greengenes and SILVA have 7 levels/columns. If you want to use all levels, adapt the Rhea script, "taxonomic_binning.R" to 7 levels and unite one more column.
             ]%>%unite(taxonomy,c("Domain","Phylum","Class","Order","Family","Genus"), sep=";",remove=TRUE) #<-Adjust accordingly if "species" taxonomic levels should be included!
  R1$taxonomy<-paste0("",R1$taxonomy,sep=";")
  R2 <- sapply(R1, gsub, pattern = "NA", replacement = "", fixed = TRUE)
  R3<-cbind(otu_table(physeq),R2) #<- Or replace wih with the name of your phyloseq object!
  R4 <- gsub(";","/t", R3[,1:16])#<-Number of samples counting from 1 + 1
  R4 <- cbind(rownames(R3), data.frame(R3, row.names=NULL))
  names(R4)[1]<-"#OTUId"
  write.table(R4,"OTU_table_w_taxonomy.txt",row.names=FALSE, sep= "\t") #<- Change the name of the output file if you like!
  message("Writing Rhea_OTU_Table in the file ?OTU_table_w_taxonomy.txt of your woking directory. 
          Also load your tree and your map file in this directory. Then you can start the Rhea pipeline!")
}

#Make an OTU table with taxonomy out of your physeq object
OTU_TAX_Rhea <- cbind.data.frame(otu_table(physeq1),tax_table(physeq1)) #<- Replace GlobalPatterns with the name of your phyloseq object!
#Write the OTU table with taxonomy suitable for Rhea
Physeq_to_Rhea(OTU_TAX_Rhea)
