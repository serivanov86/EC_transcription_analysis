# the pipeline requires several R packages from CRAN and Bioconductor
#loading the packages

library(data.table)
library(illuminaHumanv4.db)
library(illuminaHumanv3.db)
library(hgu133a.db)
library(hgu133acdf)
library(genefilter)
library(affy)
library(limma)
library(pvclust)
library(RColorBrewer)
library(gplots)
library(org.Hs.eg.db)
library(annotate)
library(gage)
library(enrichR)

setwd("C:/Working_directory") # set working directory

# identification of KEGG pathways, differentially regulated between groups of elite controllers and healthy controls
# the corresponding code is written for the comparison of pathways between EC group 1 and healthy controls
# to calculate differentially regulated pathways for other groups of elite controllers and cART-treated progressors, replace "EC1" in the code to "EC2", "EC3", "EC4", "EC5", and "HAART"


# loading data on gene symbol - KEGG pathway associations. The corresponding file can be obtained from the "Enrichr" resource (https://maayanlab.cloud/Enrichr)

gs<-readList("./Pathway data/KEGG_2019_Human")

# converting of gene symbols to Entrez IDs

data(egSymb)
gs<-sapply(gs,sym2eg)
nm<-names(gs)
gs<-sapply(1:length(gs),function(x) {gs[[x]][!is.na(gs[[x]])]})
names(gs)<-nm

# loading the table with data on 51 elite controllers, 32 cART-treated progressors, and ten healthy individuals containing only the top 50% genes with the highest variance

gse1<-fread("./Normalized data/GSE87620_normalized_filt.csv",data.table=F)
rownames(gse1)<-gse1[,1]
gse1<-gse1[,-1]

# loading of the table with information on patients' groups

groups<-read.delim("./Results/Clusters.txt",as.is=T)

# creation of the same order of sample names in tables with transcription data and data on clusters

groups<-groups[match(colnames(gse1),groups$sample),]

# creation of variables with indexes of columns from the table belonging to the EC1 group and healthy controls

groups<-groups$cluster
groups<-as.factor(groups)
group1<-grep("EC1",groups)
group2<-grep("HIVneg",groups)

# identification of KEGG pathways using gene set enrichment analysis implemented in "gage" R package

pathways <- gage(gse1, gsets = gs,
                 ref = group2, samp = group1,compare='unpaired',
                 same.dir=T)

# selection of pathways with adjusted p-value less than 0.1

pathways<-sigGeneSet(pathways,heatmap=F)

# creation of lists of up- and down-regulated pathways

pathways.up<-pathways$greater
pathways.down<-pathways$less
write.table(pathways.up,file="./Results/KEGG up EC1 - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)
write.table(pathways.down,file="./Results/KEGG down EC1 - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)

###############################################################################################################################################################################################

# identification of KEGG pathways, differentially regulated between untreated progressors in acute/chronic phase and healthy controls (GSE6740 dataset)

# loading the normalized data from the GSE6740 dataset

gse2<-fread("./Normalized data/GSE6740_normalized.csv",data.table=F)
rownames(gse2)<-gse2[,1]
gse2<-gse2[,-1]

# creation of two tables: untreated progressors in acute phase plus healthy controls ("AI1" variable), and untreated progressors in chronic phase plus healthy controls ("CI" variable)

AI1<-c("GSM155179.CEL", "GSM155181.CEL", "GSM155183.CEL", "GSM155185.CEL", "GSM155187.CEL", "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL")
CI<-c("GSM155190.CEL", "GSM155195.CEL", "GSM155201.CEL", "GSM155203.CEL", "GSM155206.CEL", "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL")
AI1<-gse2[,AI1]
CI<-gse2[,CI]

# identification of KEGG pathways associated with untreated progressors in acute phase using gene set enrichment analysis implemented in "gage" R package

pathways <- gage(AI1, gsets = gs,
                 ref = 6:10, samp = 1:5,compare='unpaired',
                 same.dir=T)

# selection of pathways with adjusted p-value less than 0.1

pathways<-sigGeneSet(pathways,heatmap=F)

# creation of lists of up- and down-regulated pathways associated with untreated progressors in the acute phase

pathways.up<-pathways$greater
pathways.down<-pathways$less
write.table(pathways.up,file="./Results/KEGG up AI1 - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)
write.table(pathways.down,file="./Results/KEGG down AI1 - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)

# identification of KEGG pathways associated with untreated progressors in  chronic phase using gene set enrichment analysis implemented in "gage" R package

pathways <- gage(CI, gsets = gs,
                 ref = 6:10, samp = 1:5,compare='unpaired',
                 same.dir=T)

# selection of pathways with adjusted p-value less than 0.1

pathways<-sigGeneSet(pathways,heatmap=F)

# creation of lists of up- and down-regulated pathways associated with untreated progressors in chronic phase

pathways.up<-pathways$greater
pathways.down<-pathways$less
write.table(pathways.up,file="./Results/KEGG up CI - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)
write.table(pathways.down,file="./Results/KEGG down CI - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)

###############################################################################################################################################################################################


# identification of KEGG pathways, differentially regulated between untreated progressors in the acute phase and healthy controls (GSE25669 dataset)

# loading the normalized data from the GSE25669 dataset

gse3<-fread("./Normalized data/GSE25669_normalized.csv",data.table=F)
rownames(gse3)<-gse3[,1]
gse3<-gse3[,-1]

# identification of KEGG pathways associated with untreated progressors in acute phase using gene set enrichment analysis implemented in "gage" R package

pathways <- gage(gse3, gsets = gs,
                 ref = 5:6, samp = 1:4,compare='unpaired',
                 same.dir=T)

# selection of pathways with adjusted p-value less than 0.1

pathways<-sigGeneSet(pathways,heatmap=F)

# creation of lists of up- and down-regulated pathways associated with untreated progressors in the acute phase

pathways.up<-pathways$greater
pathways.down<-pathways$less
write.table(pathways.up,file="./Results/KEGG up AI2 - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)
write.table(pathways.down,file="./Results/KEGG down AI2 - HIVneg.txt", row.names=T,sep='\t',quote=FALSE)
