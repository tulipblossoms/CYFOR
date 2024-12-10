library(ggplot2)
library(dplyr)
library(bioMart)
library(pheatmap)

#reading in files
myc_wt_oscillating<-read.csv("C:/Users/kelly/Internship/clock_project/ML/data/infile_Shep_MYC_OFF.csv")
myc_oe_oscillating<-read.csv("C:/Users/kelly/Internship/clock_project/ML/data/infile_Shep_MYC_ON.csv")
myc_combined<-read.csv("C:/Users/kelly/Internship/clock_project/ML/data/myc_combined_list_filtered_v2.csv")
myc_wt_oscillating[,2:ncol(myc_wt_oscillating)]<-apply(myc_wt_oscillating[,2:ncol(myc_wt_oscillating)],2,as.numeric)
myc_oe_oscillating[,2:ncol(myc_oe_oscillating)]<-apply(myc_oe_oscillating[,2:ncol(myc_oe_oscillating)],2,as.numeric)

myc_combined<-cbind(myc_wt_oscillating,myc_oe_oscillating[,2:ncol(myc_oe_oscillating)])
#myc_wt_oscillating<-log2(myc_wt_oscillating+1)
rownames(myc_combined)<-myc_combined$ENSG
myc_combined<-myc_combined[,2:ncol(myc_combined)]

#PERFORM QUANTILE NORMALIZATION
myc_normalized<-qn(myc_combined)
rownames(myc_wt_oscillating_normalzed)<-rownames(myc_wt_oscillating)
colnames(myc_wt_oscillating_normalzed)<-colnames(myc_wt_oscillating)
oscillating_genes<-read.csv("C:/Users/kelly/Internship/clock_project/ML/data/SHEPMYCOFFOscillatingGenes.csv")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(myc_wt_oscillating),
  mart = mart
)

annotations <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = oscillating_genes$x,
  mart = mart
)

annotation1<-which(annotations$ensembl_gene_id%in%myc_wt_oscillating$ENSG)
ensembl_id<-annotations$ensembl_gene_id[annotation1]
myc_wt_oscillating_normalzed<-myc_combined[which(rownames(myc_combined)%in%ensembl_id),]
myc_wt_oscillating_normalzed<-log2(myc_wt_oscillating_normalzed+1)
normalized_off<-myc_wt_oscillating_normalzed[,1:28]
normalized_on<-myc_wt_oscillating_normalzed[,29:56]
normalized_off<-apply(normalized_off,2,as.numeric)
normalized_on<-apply(normalized_on,2,as.numeric)


names<-rownames(myc_combined_off)
myc_combined_off<-apply(myc_combined_off,2,as.numeric)
for(i in 1:nrow(myc_combined_off)){
  row_mean<-mean(myc_combined_off[i,])
  row_sd<-sd(myc_combined_off[i,])
  myc_combined_off[i,]<-(myc_combined_off[i,]-row_mean)/row_sd
}
names<-rownames(myc_combined_on)
myc_combined_on<-apply(myc_combined_on,2,as.numeric)
for(i in 1:nrow(myc_combined_on)){
  row_mean<-mean(myc_combined_on[i,])
  row_sd<-sd(myc_combined_on[i,])
  myc_combined_on[i,]<-(myc_combined_on[i,]-row_mean)/row_sd
}
rownames(myc_combined_on)<-names
colors<-c("#4F97C1","#70AFD0", "#91C8DF" ,"#B2E1EE", "#D3F9FD", "#FCF1E4","#FDBE8F", "#FD9963","#E75119", "#FC7436", "#B93D13", "#8B2911")
intervals<-c(-3.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,2,4)
annotation_col <- data.frame(
  Group = c(rep("Sample X",14),rep("Sample Y",14)),
  Day=c(rep("Day 1",7),rep("Day 2",7),rep("Day 1",7),rep("Day 2",7))
)
rownames(annotation_col) <- colnames(myc_combined_on)

map<-pheatmap(myc_combined_off[,15:28],cluster_cols=FALSE,show_rownames=FALSE,fontsize=15,color=colors,breaks=intervals,main="SHEP MYC OFF Oscillating Genes",annotation_col=annotation_col,annotation_names_col=FALSE)
rownames(annotation_col) <- colnames(myc_combined_on)

map1<-pheatmap(myc_combined_on[map$tree_row$order,15:28],cluster_rows=FALSE,fontsize=15,cluster_cols=FALSE,show_rownames=FALSE,color=colors,breaks=intervals,main="SHEP MYC ON Oscillating Genes",annotation_col=annotation_col,annotation_names_col=FALSE)
