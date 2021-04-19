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

# identification of genes differentially expressed between groups of elite controllers and healthy controls
# the corresponding code is written for comparison of transcription between EC group 1 and healthy controls
# to calculate differentially expressed genes for other groups of elite controllers and cART-treated progressors, replace "EC1" in the code to "EC2", "EC3", "EC4", "EC5", and "HAART"

# loading the normalized data on 51 elite controllers, 32 cART-treated progressors, and ten healthy individuals

tt<-fread("./Normalized data/GSE87620_normalized_filt.csv",data.table=F)
rownames(tt)<-tt[,1]
tt<-tt[,-1]

# creation of the same order of sample names in tables with transcription data and data on clusters

groups<-read.delim("./Results/Clusters.txt",as.is=T)
groups<-groups[match(colnames(tt),groups$sample),]
groups<-groups$cluster

# selection of columns, which correspond to EC group 1 and healthy controls

tt<-tt[,groups %in% c("EC1","HIVneg")]

# creation of "groups" variable with data on sample group: EC group 1 or healthy control

groups<-groups[groups %in% c("EC1","HIVneg")]
groups<-as.factor(groups)

# performing the estimation of differentially expressed genes using the LIMMA method

design <- model.matrix(~0+groups,tt)
rownames(design) <- colnames(tt)
colnames(design) <- c("EC1","HIVneg")
fit <- lmFit(tt,design)
contrast.matrix <- makeContrasts(EC1-HIVneg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef="EC1 - HIVneg", n=500000, sort.by="p") # saving data on all genes (n=500000) sorted by p-value

# adding data on gene symbols to the table with differentially expressed genes based on Entrez IDs. Saving data to file

deg <- data.frame(EntrezID = rownames(deg),Symbol = getSYMBOL(rownames(deg), data='org.Hs.eg'),deg,stringsAsFactors = F)
deg<-na.omit(deg)
write.table(deg,file="./Results/EC1 - HIVneg.txt", row.names=F,col.names=T,sep='\t',quote=FALSE)

###########################################################################################################################################################


# identification of genes differentially expressed between untreated progressors in acute and chronic phases, and healthy controls (GSE6740 dataset)

# loading the normalized data from the GSE6740 dataset

tt<-fread("./Normalized data/GSE6740_normalized.csv",data.table=F)
rownames(tt)<-tt[,1]
tt<-tt[,-1]

# creation of two tables: untreated progressors in acute phase plus healthy controls ("AI1" variable), and untreated progressors in chronic phase plus healthy controls ("CI" variable)

AI1<-c("GSM155179.CEL", "GSM155181.CEL", "GSM155183.CEL", "GSM155185.CEL", "GSM155187.CEL", "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL")
CI<-c("GSM155190.CEL", "GSM155195.CEL", "GSM155201.CEL", "GSM155203.CEL", "GSM155206.CEL", "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL")
AI1<-tt[,AI1]
CI<-tt[,CI]

# performing the estimation of differentially expressed genes using the LIMMA method for the untreated progressors in acute phase

groups<-as.factor(c(rep("AI1",5),rep("HIVneg",5)))
design <- model.matrix(~0+groups,AI1)
rownames(design) <- colnames(AI1)
colnames(design) <- c("AI1","HIVneg")
fit <- lmFit(AI1,design)
contrast.matrix <- makeContrasts(AI1-HIVneg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef="AI1 - HIVneg", n=500000, sort.by="p") # saving data on all genes (n=500000) sorted by p-value

# adding data on gene symbols to the table with differentially expressed genes based on Entrez IDs. Saving data to file

deg <- data.frame(EntrezID = rownames(deg),Symbol = getSYMBOL(rownames(deg), data='org.Hs.eg'),deg,stringsAsFactors = F)
deg<-na.omit(deg)
write.table(deg,file="./Results/AI1 - HIVneg.txt", row.names=F,col.names=T,sep='\t',quote=FALSE)

# performing the estimation of differentially expressed genes using the LIMMA method for the untreated progressors in chronic phase

groups<-as.factor(c(rep("CI",5),rep("HIVneg",5)))
design <- model.matrix(~0+groups,CI)
rownames(design) <- colnames(CI)
colnames(design) <- c("CI","HIVneg")
fit <- lmFit(CI,design)
contrast.matrix <- makeContrasts(CI-HIVneg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef="CI - HIVneg", n=500000, sort.by="p") # saving data on all genes (n=500000) sorted by p-value

# adding data on gene symbols to the table with differentially expressed genes based on Entrez IDs. Saving data to file

deg <- data.frame(EntrezID = rownames(deg),Symbol = getSYMBOL(rownames(deg), data='org.Hs.eg'),deg,stringsAsFactors = F)
deg<-na.omit(deg)
write.table(deg,file="./Results/CI - HIVneg.txt", row.names=F,col.names=T,sep='\t',quote=FALSE)

###########################################################################################################################################################


# identification of genes differentially expressed between untreated progressors in the acute phase and healthy controls (GSE25669 dataset)

# loading the normalized data from the GSE25669 dataset

AI2<-fread("./Normalized data/GSE25669_normalized.csv",data.table=F)
rownames(AI2)<-AI2[,1]
AI2<-AI2[,-1]

# creation of table with transcription data on untreated progressors in acute phase plus healthy controls ("AI2" variable)

AI2<-AI2[,c("GSM630684.Signal","GSM630685.Signal","GSM630703.Signal","GSM630706.Signal",
            "GSM630701.Signal","GSM630716.Signal")]

# performing the estimation of differentially expressed genes using the LIMMA method for untreated progressors in the acute phase

groups<-as.factor(c("AI2", "AI2", "AI2","AI2","HIVneg","HIVneg"))
design <- model.matrix(~0+groups,AI2)
rownames(design) <- colnames(AI2)
colnames(design) <- c("AI2","HIVneg")
fit <- lmFit(AI2,design)
contrast.matrix <- makeContrasts(AI2-HIVneg,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef="AI2 - HIVneg", n=500000, sort.by="p") # saving data on all genes (n=500000) sorted by p-value

# adding data on gene symbols to the table with differentially expressed genes based on Entrez IDs. Saving data to file

deg <- data.frame(EntrezID = rownames(deg),Symbol = getSYMBOL(rownames(deg), data='org.Hs.eg'),deg,stringsAsFactors = F)
deg<-na.omit(deg)
write.table(deg,file="./Results/AI2 - HIVneg.txt", row.names=F,col.names=T,sep='\t',quote=FALSE)
