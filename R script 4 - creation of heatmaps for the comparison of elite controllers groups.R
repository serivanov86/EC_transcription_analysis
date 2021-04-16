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

# loading the data on differentially expressed genes for each group of elite controllers, cART-treated and untreated progressors

EC1<-read.delim("./Results/EC1 - HIVneg.txt",as.is=T)
EC2<-read.delim("./Results/EC2 - HIVneg.txt",as.is=T)
EC3<-read.delim("./Results/EC3 - HIVneg.txt",as.is=T)
EC4<-read.delim("./Results/EC4 - HIVneg.txt",as.is=T)
EC5<-read.delim("./Results/EC5 - HIVneg.txt",as.is=T)
HAART<-read.delim("./Results/HAART - HIVneg.txt",as.is=T)
AI1<-read.delim("./Results/AI1 - HIVneg.txt",as.is=T)
CI<-read.delim("./Results/CI - HIVneg.txt",as.is=T)
AI2<-read.delim("./Results/AI2 - HIVneg.txt",as.is=T)

# selection of the genes with log fold change more than |0.7| and adjusted p-value less than 0.1

EC1<-EC1[abs(EC1$logFC)>0.7 & EC1$adj.P.Val<0.1,]
EC2<-EC2[abs(EC2$logFC)>0.7 & EC2$adj.P.Val<0.1,]
EC3<-EC3[abs(EC3$logFC)>0.7 & EC3$adj.P.Val<0.1,]
EC4<-EC4[abs(EC4$logFC)>0.7 & EC4$adj.P.Val<0.1,]
EC5<-EC5[abs(EC5$logFC)>0.7 & EC5$adj.P.Val<0.1,]
HAART<-HAART[abs(HAART$logFC)>0.7 & HAART$adj.P.Val<0.1,]
AI1<-AI1[abs(AI1$logFC)>0.7 & AI1$adj.P.Val<0.1,]
CI<-CI[abs(CI$logFC)>0.7 & CI$adj.P.Val<0.1,]
AI2<-AI2[abs(AI2$logFC)>0.7 & AI2$adj.P.Val<0.1,]

# creation of "deg" variable containing a list of Entrez IDs for genes that are differentially expressed in at least one EC group or cART group compared to the healthy control 

deg<-unique(c(EC1$EntrezID,EC2$EntrezID,EC3$EntrezID,EC4$EntrezID,EC5$EntrezID,HAART$EntrezID))

# loading the table with data on 51 elite controllers, 32 cART-treated progressors, and ten healthy individuals containing only the top 50% genes with the highest variance

tt<-fread("./Normalized data/GSE87620_normalized_filt.csv",data.table=F)
rownames(tt)<-tt[,1]
tt<-tt[,-1]

# rewriting the table "tt". Saving only differentially expressed genes

tt<-tt[rownames(tt) %in% deg,]

# loading of the table with information on patients' groups

clusters<-read.delim("./Results/Clusters.txt",as.is=T)

# creation of the same order of sample names in tables with transcription data and data on clusters

clusters<-clusters[match(colnames(tt),clusters$sample),]

# creation of the palette for the heatmap

rbpal<-colorRampPalette(brewer.pal(10, "RdBu"))(256)
rbpal<-rev(rbpal)

# creation of function for the clustering of samples based on the Ward method

hclustfun<-function(x) hclust(x, method="ward.D2")

# creation of function for the calculation of distances between samples based on 1 - Pearson correlation coefficient values

distfun<-function(x) as.dist((1-cor(t(x)))/2)

# creation of the heatmap to compare the five EC groups with cART-treated progressors and healthy controls (Figure 2 in the Manuscript)

heatmap.2(data.matrix(tt), col=rbpal, trace="none",density.info="none",
          ColSideColors=as.character(clusters$color),labRow="",labCol="",key.title=NA,
          margins = c(0.5, 0.5),
          hclustfun=hclustfun,
          distfun=distfun,scale="row")

# creation of "deg" variable containing a list of Entrez IDs for genes that are differentially expressed in at least one EC group, cART-treated progressors, untreated progressors in acute (AI1, AI2) and chronic (CI) phases compared to the healthy control

deg<-unique(c(EC1$EntrezID,EC2$EntrezID,EC3$EntrezID,EC4$EntrezID,EC5$EntrezID,HAART$EntrezID,AI1$EntrezID,CI$EntrezID,AI2$EntrezID))

# loading the corresponding normalized transcription data (GSE87620, GSE6740 and GSE25669 datasets)

gse1<-fread("./Normalized data/GSE87620_normalized_filt.csv",data.table=F)
rownames(gse1)<-gse1[,1]
gse1<-gse1[,-1]

gse2<-fread("./Normalized data/GSE6740_normalized.csv",data.table=F)
rownames(gse2)<-gse2[,1]
gse2<-gse2[,-1]

gse3<-fread("./Normalized data/GSE25669_normalized.csv",data.table=F)
rownames(gse3)<-gse3[,1]
gse3<-gse3[,-1]

# creation of "deg" variable containing a list of Entrez IDs for genes that are differentially expressed in at least one of investigated groups

deg<-Reduce(intersect, list(deg,rownames(gse1),rownames(gse2),rownames(gse3)))

# loading of the table with information on patients' groups

clusters<-read.delim("./Results/Clusters.txt",as.is=T)

# creation of the same order of sample names in tables with transcription data on elite controllers and cART-treated progressors ("gse1" variable) and data on groups of patients

clusters<-clusters[match(colnames(gse1),clusters$sample),]

# creation of variables with log-fold changes for each EC group and cART-treated progressors

EC1.fold<-apply(gse1[,clusters$cluster=="EC1"],1,mean)-apply(gse1[,clusters$cluster=="HIVneg"],1,mean)
EC2.fold<-apply(gse1[,clusters$cluster=="EC2"],1,mean)-apply(gse1[,clusters$cluster=="HIVneg"],1,mean)
EC3.fold<-apply(gse1[,clusters$cluster=="EC3"],1,mean)-apply(gse1[,clusters$cluster=="HIVneg"],1,mean)
EC4.fold<-apply(gse1[,clusters$cluster=="EC4"],1,mean)-apply(gse1[,clusters$cluster=="HIVneg"],1,mean)
EC5.fold<-apply(gse1[,clusters$cluster=="EC5"],1,mean)-apply(gse1[,clusters$cluster=="HIVneg"],1,mean)
HAART.fold<-apply(gse1[,clusters$cluster=="HAART"],1,mean)-apply(gse1[,clusters$cluster=="HIVneg"],1,mean)

# creation of two separate tables with transcription data from untreated progressors in acute (AI1) and chronic (CI) phases (GSE6740 dataset)

AI1<-c("GSM155179.CEL", "GSM155181.CEL", "GSM155183.CEL", "GSM155185.CEL", "GSM155187.CEL", "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL")
CI<-c("GSM155190.CEL", "GSM155195.CEL", "GSM155201.CEL", "GSM155203.CEL", "GSM155206.CEL", "GSM155229.CEL", "GSM155232.CEL", "GSM155234.CEL", "GSM155236.CEL", "GSM155238.CEL")
AI1<-gse2[,AI1]
CI<-gse2[,CI]

# creation of variables with log-fold changes for untreated progressors in acute (AI1.fold) and chronic (CI.fold) phases (GSE6740 dataset)

groups<-as.factor(c(rep("AI1",5),rep("HIVneg",5)))
AI1.fold<-apply(AI1[,groups=="AI1"],1,mean)-apply(AI1[,groups=="HIVneg"],1,mean)

groups<-as.factor(c(rep("CI",5),rep("HIVneg",5)))
CI.fold<-apply(CI[,groups=="CI"],1,mean)-apply(CI[,groups=="HIVneg"],1,mean)

# creation of table with transcription data from untreated progressors in the acute phase (AI2) (GSE25669 dataset)

AI2<-gse3[,c("GSM630684.Signal","GSM630685.Signal","GSM630703.Signal","GSM630706.Signal",
             "GSM630701.Signal","GSM630716.Signal")]

# creation of variable with log-fold changes for untreated progressors in the acute phase (AI2.fold) (GSE25669 dataset)

groups<-as.factor(c("AI2", "AI2", "AI2","AI2","HIVneg","HIVneg"))
AI2.fold<-apply(AI2[,groups=="AI2"],1,mean)-apply(AI2[,groups=="HIVneg"],1,mean)

# rewriting the variables with log-fold changes. The new variables contain only genes that are differentially expressed in at least one of the patients' groups

EC1.fold<-EC1.fold[names(EC1.fold) %in% deg]
EC2.fold<-EC2.fold[names(EC2.fold) %in% deg]
EC3.fold<-EC3.fold[names(EC3.fold) %in% deg]
EC4.fold<-EC4.fold[names(EC4.fold) %in% deg]
EC5.fold<-EC5.fold[names(EC5.fold) %in% deg]
HAART.fold<-HAART.fold[names(HAART.fold) %in% deg]
AI1.fold<-AI1.fold[names(AI1.fold) %in% deg]
CI.fold<-CI.fold[names(CI.fold) %in% deg]
AI2.fold<-AI2.fold[names(AI2.fold) %in% deg]

# creation of table with log-fold changes

fold<-as.data.frame(cbind(EC1.fold,EC2.fold,EC3.fold,EC4.fold,EC5.fold,HAART.fold,AI1.fold,CI.fold,AI2.fold))
fwrite(fold,file="./Results/fold.csv",row.names=T)

# creation of the palette for the heatmap

rbpal<-colorRampPalette(brewer.pal(10, "RdBu"))(256)
rbpal<-rev(rbpal)

# creation of function for the clustering of samples based on the Ward method

hclustfun<-function(x) hclust(x, method="ward.D2")

# creation of function for the calculation of distances between samples based on 1 - Pearson correlation coefficient values

distfun<-function(x) as.dist((1-cor(t(x)))/2)

# creation of vectors with colors for heatmap's horizontal sidebar and labels for columns

pat.col<-as.character(ncol(fold))
pat.col<-c("blue","cyan","green","red","yellow","grey","black","black","black")
lb<-c("EC1","EC2","EC3","EC4","EC5","cART","AI1","CI","AI2")

# creation of the heatmap to compare fold changes of five EC groups with cART-treated progressors, and untreated progressors in acute and chronic phases and healthy controls (Figure 2 in the Manuscript)

heatmap.2(data.matrix(fold), col=rbpal, trace="none",density.info="none",
          ColSideColors=pat.col,labRow="",labCol=lb,key.title=NA,
          margins = c(4, 0.5),
          hclustfun=hclustfun,
          distfun=distfun,scale="row")
