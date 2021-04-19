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

# the corresponding code is written for the identification of pathways associated with genes, which are differentially expressed in EC group 2 but are not differentially expressed in EC group 3
# to identify pathways associated with genes, which are differentially expressed in EC group 3 but are not differentially expressed in EC group 4, replace "EC2" and "EC3" with "EC3" and "EC4"

# loading the data on differentially expressed genes for EC groups 2 and 3

EC2<-read.delim("./Results/EC2 - HIVneg.txt", as.is=T)
EC3<-read.delim("./Results/EC3 - HIVneg.txt", as.is=T)

# selection the up- and down-regulated genes with log fold change more than |0.7| and adjusted p-value less than 0.1
# the corresponding thresholds were chosen empirically to balance the number of DEGs and statistical significance of differential expression

EC2.up<-EC2$Symbol[EC2$logFC>0.7 & EC2$adj.P.Val<0.1]
EC2.down<-EC2$Symbol[EC2$logFC<(-0.7) & EC2$adj.P.Val<0.1]
EC3.up<-EC3$Symbol[EC3$logFC>0.7 & EC3$adj.P.Val<0.1]
EC3.down<-EC3$Symbol[EC3$logFC<(-0.7) & EC3$adj.P.Val<0.1]

# identification of up- and down-regulated genes, which are differentially expressed in EC group 2 but are not differentially expressed in EC group 3

DEG.up<-EC2.up[!EC2.up %in% EC3.up]
DEG.down<-EC2.down[!EC2.down %in% EC3.down]

# performing the over-representation analysis for up-regulated genes using the "enrichR" R package

enriched <- enrichr(DEG.up, "KEGG_2019_Human")

# parsing of obtained results and creation of table with KEGG pathways

DEP.up<-enriched$KEGG_2019_Human
yy<-strsplit(DEP.up$Overlap,split="/",fixed=T)
yy<-t(sapply(1:length(yy),function(x){
  yy[[x]]
}))
yy<-apply(yy,2,as.numeric)
DEP.up<-data.frame(Term=DEP.up[,1],n=yy[,1],N=yy[,2],DEP.up[,c(3:4,7:9)],stringsAsFactors = F)

# selection of pathways with adjusted p-value less than 0.1 and number of studied genes in a pathway more than 2

DEP.up<-DEP.up[DEP.up$Adjusted.P.value<0.1 & DEP.up$n>2,]
write.table(DEP.up,file="./Results/KEGG EC2-EC3 up.txt", row.names=F,sep='\t',quote=FALSE)

# performing the over-representation analysis for down-regulated genes using the "enrichR" R package

enriched <- enrichr(DEG.down, "KEGG_2019_Human")

# parsing of obtained results and creation of table with KEGG pathways

DEP.down<-enriched$KEGG_2019_Human
yy<-strsplit(DEP.down$Overlap,split="/",fixed=T)
yy<-t(sapply(1:length(yy),function(x){
  yy[[x]]
}))
yy<-apply(yy,2,as.numeric)
DEP.down<-data.frame(Term=DEP.down[,1],n=yy[,1],N=yy[,2],DEP.down[,c(3:4,7:9)],stringsAsFactors = F)

# selection of pathways with adjusted p-value less than 0.1 and number of studied genes in a pathway more than 2

DEP.down<-DEP.down[DEP.down$Adjusted.P.value<0.1 & DEP.down$n>2,]
write.table(DEP.down,file="./Results/KEGG EC2-EC3 down.txt", row.names=F,sep='\t',quote=FALSE)
