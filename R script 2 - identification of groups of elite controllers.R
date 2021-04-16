
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
dir.create("./Results") # create a folder for results of pipeline at intermediate steps

# loading the normalized data on 51 elite controllers, 32 cART-treated progressors, and ten healthy individuals

tt<-fread("./Normalized data/GSE87620_normalized.csv",data.table=F)
rownames(tt)<-tt[,1]
tt<-tt[,-1]

# selection of columns of "tt" table, which correspond to 51 elite controllers. Rewritting the variable "tt"

inx <- grep("EC",colnames(tt))
tt<-tt[,inx]

# identification of Entrez IDs ("nm" variable) for the top 50% of genes with the highest variance

vr<-apply(tt,1,var)
vr<-sort(vr,decreasing=T)
nm<-names(vr)[1:floor(length(vr)/2)]

# rewriting the table with normalized transcription data. Saving only the top 50% of genes with the highest variance

tt <- tt[nm,]
rownames(tt)<-nm
fwrite(tt,file="./Results/EC_filt.csv",row.names=T)

# creation of function for the clustering of samples based on the Ward method

hclustfun<-function(x) hclust(x, method="ward.D2")

# creation of function for the calculation of distances between samples based on 1 - Pearson correlation coefficient values

distfun<-function(x) as.dist((1-cor(t(x)))/2)

# cluster analysis and bootstrap resampling analysis for the estimation of clusters' stability

fit <- pvclust(tt, method.hclust="ward.D2",parallel=T,nboot=1000)
plot(fit) # creation of dendrogram

cl<-fit$hclust
clusters<-cutree(cl,k=5) # selection of five clusters
table(clusters) # printing of clusters' sizes

# creation of table with information on clusters

clusters<-data.frame(sample=colnames(tt),cluster=clusters,color=character(length(clusters)),stringsAsFactors = F)
clusters<-clusters[order(clusters$cluster),]
clusters$cluster<-as.character(clusters$cluster)
clusters$cluster[1:5]<-"EC2"
clusters$cluster[6:21]<-"EC5"
clusters$cluster[22:30]<-"EC3"
clusters$cluster[31:38]<-"EC1"
clusters$cluster[39:51]<-"EC4"
clusters$color[1:5]<-"cyan"
clusters$color[6:21]<-"yellow"
clusters$color[22:30]<-"green"
clusters$color[31:38]<-"blue"
clusters$color[39:51]<-"red"
write.table(clusters,file="./Results/Clusters.txt", row.names=F,col.names=T,sep='\t',quote=FALSE)


# creation of the same order of sample names in tables with transcription data and data on clusters

clusters<-clusters[match(colnames(tt),clusters$sample),]

# creation of the palette for the heatmap

rbpal <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
rbpal <- rev(rbpal)

# creation of the heatmap to visualize the clusters (Figure 1 in the Manuscript)

heatmap.2(data.matrix(tt), col=rbpal, trace="none",density.info="none",
          ColSideColors=as.character(clusters$color),labRow="",labCol="",key.title=NA,
          margins = c(0.5, 0.5),
          hclustfun=hclustfun,
          distfun=distfun,scale="row")

# loading the normalized data on 51 elite controllers, 32 cART-treated progressors, and ten healthy individuals

tt<-fread("./Normalized data/GSE87620_normalized.csv",data.table=F)
rownames(tt)<-tt[,1]
tt<-tt[,-1]

# loading the normalized transcription data on 51 elite controllers with top 50% genes with the highest variance

tt.EC<-fread("./Results/EC_filt.csv",data.table=F)
rownames(tt.EC)<-tt.EC[,1]
tt.EC<-tt.EC[,-1]

# rewriting the table with normalized transcription data. Saving only the top 50% of genes with the highest variance

nm<-rownames(tt.EC)
tt<-tt[nm,]
rownames(tt)<-nm
fwrite(tt,file="./Normalized data/GSE87620_normalized_filt.csv",row.names=T)

# adding the data on cART-treated patients and healthy people to the "Clusters" table

nm<-colnames(tt)
clusters<-read.delim("./Results/Clusters.txt",as.is=T)
inx <- c(grep("HAART",nm),grep("HIVneg",nm))
clusters.add<-data.frame(sample=nm[inx],cluster=character(length(nm[inx])),color=character(length(nm[inx])),stringsAsFactors = F)
clusters.add$cluster[1:32]<-"HAART"
clusters.add$cluster[33:42]<-"HIVneg"
clusters.add$color[1:32]<-"grey"
clusters.add$color[33:42]<-"black"
clusters<-rbind(clusters,clusters.add)
write.table(clusters,file="./Results/Clusters.txt", row.names=F,col.names=T,sep='\t',quote=FALSE)
