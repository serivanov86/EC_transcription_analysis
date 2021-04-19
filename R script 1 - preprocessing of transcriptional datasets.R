
# the pipeline requires several R packages from CRAN and Bioconductor

# installation of packages

install.packages("data.table")
install.packages("pvclust")
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("enrichR")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("illuminaHumanv3.db")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133acdf")
BiocManager::install("genefilter")
BiocManager::install("limma")
BiocManager::install("affy")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("annotate")
BiocManager::install("gage")

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
dir.create("./Normalized data") # create a folder for normalized data

# preprocessing of the data on 51 elite controllers, 32 cART-treated progressors, and 10 healthy individuals
# GSE87620 dataset was created by Illumina HumanHT-12 V4.0 platform

# loading the data from the file, which was obtained from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87620)

tt<-fread("./Unnormalized data from GEO/GSE87620_non_normalized.txt",data.table=F)
rownames(tt)<-tt$ID_REF
tt<-tt[,-1]

# creation of EListRaw object, which is a list-based S4 class for storing expression values

nm<-colnames(tt)
j<-seq(1,length(nm)-1,by=2)
unnorm.data <- new('EListRaw')
unnorm.data$E<-data.matrix(tt[,j])
unnorm.data$other$Detection<-data.matrix(tt[,j+1])

# normalization of data

norm.data<-neqc(unnorm.data)

# removing low-expressed probes

expressed <- rowSums(norm.data$other$Detection < 0.05) >= 1
norm.data <- norm.data[expressed,]
eset<-norm.data$E

# removing probes without Entrez ID and selecting one probe per gene with the largest variance across samples

eset<-ExpressionSet(assayData=eset,annotation="illuminaHumanv4")
eset.filt<-nsFilter(eset,require.entrez=T, remove.dupEntrez=T, var.filter=F)
eset.filt<-eset.filt$eset
eset.filt<-exprs(eset.filt)

# changing probe IDs to Entrez IDs and saving the normalized matrix to file

EntrezID <- unlist(mget(rownames(eset.filt), illuminaHumanv4ENTREZID))
rownames(eset.filt)<-EntrezID
fwrite(as.data.frame(eset.filt),file="./Normalized data/GSE87620_normalized.csv",row.names=T)

###############################################################################

# preprocessing of the data on four untreated progressors in the acute phase, and two healthy individuals
# GSE25669 dataset was created by Illumina HumanHT-12 V3.0 platform

# loading the data from the file, which was obtained from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25669)

tt<-fread("./Unnormalized data from GEO/GSE25669_non_normalized.txt",data.table=F)
rownames(tt)<-tt$ID_REF
tt<-tt[,-1]

# creation of EListRaw object, which is a list-based S4 class for storing expression values

nm<-colnames(tt)
j<-seq(1,length(nm)-1,by=2)
unnorm.data <- new('EListRaw')
unnorm.data$E<-data.matrix(tt[,j])
unnorm.data$other$Detection<-data.matrix(tt[,j+1])

# normalization of data

norm.data<-neqc(unnorm.data)

# removing low-expressed probes

expressed <- rowSums(norm.data$other$Detection < 0.05) >= 1
norm.data <- norm.data[expressed,]
eset<-norm.data$E

# removing probes without Entrez ID and selecting one probe per gene with the largest variance across samples

eset<-ExpressionSet(assayData=eset,annotation="illuminaHumanv3")
eset.filt<-nsFilter(eset,require.entrez=T, remove.dupEntrez=T, var.filter=F)
eset.filt<-eset.filt$eset
eset.filt<-exprs(eset.filt)

# changing probeIDs to EntrezIDs, and saving the normalized matrix to file

EntrezID <- unlist(mget(rownames(eset.filt), illuminaHumanv3ENTREZID))
rownames(eset.filt)<-EntrezID
fwrite(as.data.frame(eset.filt),file="./Normalized data/GSE25669_normalized.csv",row.names=T)

###############################################################################

# preprocessing of the data on five untreated progressors in the acute phase, five untreated progressors in the chronic phase, and five healthy individuals
# GSE6740 dataset was created by Affymetrix Human Genome U133A platform
# loading the data from CEL files, which were obtained from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6740)
# Files for acute HIV infection: "GSM155179.CEL", "GSM155181.CEL", "GSM155183.CEL", "GSM155185.CEL", "GSM155187.CEL"
# Files for chronic HIV infection: "GSM155190.CEL", "GSM155195.CEL", "GSM155201.CEL", "GSM155203.CEL", "GSM155206.CEL"
# Files for healthy people: "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL"
# These files must be located in the working directory

# loading the data from CEL files

cels <- list.files("./Unnormalized data from GEO/", pattern="CEL")
tt <- ReadAffy(filenames = paste("./Unnormalized data from GEO/",cels,sep=""), cdfname='hgu133acdf')

# normalization of data

eset <- rma(tt)
eset <- exprs(eset)

# removing low-expressed probes

calls<-mas5calls(tt)
calls<-exprs(calls)
nm<-colnames(calls)
calls<-sapply(1:ncol(calls),function(x) {ifelse(calls[,x]=="A",1,0)})
colnames(calls)<-nm
sm<-apply(calls,1,sum)
eset<-eset[sm!=ncol(eset),]

# removing probes without Entrez ID and selecting one probe per gene with the largest variance across samples

eset<-ExpressionSet(assayData=eset,annotation="hgu133a")
eset.filt<-nsFilter(eset,require.entrez=T, remove.dupEntrez=T, var.filter=F)
eset.filt<-eset.filt$eset
eset.filt<-exprs(eset.filt)

# changing probe IDs to Entrez IDs and saving the normalized matrix to file

EntrezID <- unlist(mget(rownames(eset.filt), hgu133aENTREZID))
rownames(eset.filt)<-EntrezID
fwrite(as.data.frame(eset.filt),file="./Normalized data/GSE6740_normalized.csv",row.names=T)
