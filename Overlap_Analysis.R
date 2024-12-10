library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)
library(dplyr)
library(R.utils)

#DELETE THE URL AT THE TOP
metadata<-read.table("C:/Users/kelly/Internship/clock_project/ML/data/experiment_report_2024_11_5_23h_57m.tsv",header=TRUE,sep="\t",fill=TRUE,quote="")
metadata<-metadata[which(metadata$Organism=="Homo sapiens"),]
metadata<-metadata[which(metadata$Target.gene.symbol%in%gene_metadata[,2]),]
metadata_unique <- metadata %>%
  distinct(Target.gene.symbol, .keep_all = TRUE)


files<-readLines("C:/Users/kelly/Internship/clock_project/ML/data/files.txt")
files<-files[grepl("\\.bed\\.gz$",files)]

regulators<-matrix(0,1004,978)
symbols<-gene_metadata[which(gene_metadata[,1]%in%rownames(regulators)),2]
rownames(regulators)<-unique(substr(combined_list$X,1,18))
colnames(regulators)<-metadata_unique$Target.of.assay
for(j in 1:length(files)){
  target_file<-substr(files[j],nchar(files[j])-17,nchar(files[j]))
  target_tf_name<-metadata_unique$Target.gene.symbol[which(grepl(substr(files[j],nchar(files[j])-17,nchar(files[j])-7),metadata_unique$Files))]
  if(length(target_tf_name)!=0){
    target_tf<-gene_metadata[which(gene_metadata[,2]==target_tf_name),1]
    bed_data <- rtracklayer::import(paste0("C:/Users/kelly/Internship/clock_project/ML/data/tf_data/kelly_tf/",target_file), 
                                    format = "BED", 
                                    extraCols = c(
                                      signalValue = "numeric",
                                      pValue = "numeric",
                                      qValue = "numeric",
                                      peak = "integer"
                                    ))
    bed_data<-bed_data[order(bed_data$qValue,decreasing=FALSE),]
    bed_data_df<-as.data.frame(bed_data)


    peak_annotations <- annotatePeak(bed_data, 
                                     TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                     annoDb = "org.Hs.eg.db")
    regulators[,which(colnames(regulators)==target_tf_name)]<-ifelse(rownames(regulators)%in%unique(peak_annotations@anno$SYMBOL),1,0)
}
}

genes<-unique(substr(rownames(regulators),1,18))
gene_prefixes <- substr(rownames(regulators), 1, 18)

new_matrix<-matrix(0,54464,2928)
for (j in seq_along(colnames(regulators))) {
  non_zero_genes <- genes[regulators[gene_indices, j] != 0]
  new_regulators_ind <- seq(j, j + 2)
  for (gene in non_zero_genes) {
    ind <- which(substr(combined_list$X, 1, 18) == gene)
    replacement_matrix <- combined_list[ind, c(2, 7, 14)]
    replacement_matrix[,3]<-abs(replacement_matrix[,3])
    new_matrix[which(substr(rownames(new_matrix),1,18)==gene),new_regulators_ind]<-apply(replacement_matrix,2,as.numeric)
  }
}

new_matrix<-rep(0,nrow(ML_dataset))
for(tf in 1:length(unique_tfs)){
  A_indiv<-rep(0,nrow(ML_dataset))
  B0_indiv<-rep(0,nrow(ML_dataset))
  phase_indiv<-rep(0,nrow(ML_dataset))
  non_zero_genes <- genes[regulators[gene_indices, unique_tfs[tf]] != 0]
  non_zero_genes_indexes <- regulators[gene_indices, unique_tfs[tf]] != 0
  A_values<-combined_list_train$A_val[ordered[tf]:(ordered[tf]+63)]
  B0_values<-combined_list_train$B0[ordered[tf]:(ordered[tf]+63)]
  phase_values<-combined_list_train$hour[ordered[tf]:(ordered[tf]+63)]
  if(length(non_zero_genes)>0){
    for (gene in 1:length(non_zero_genes)) {
      if (non_zero_genes_indexes[gene]==TRUE) {  # Only proceed if there are matches
        A_indiv[gene_index_list[[gene]]]<-A_values
        B0_indiv[gene_index_list[[gene]]] <- B0_values
        phase_indiv[gene_index_list[[gene]]] <- abs(phase_values)
      }
    }
  }
  
  new_matrix<-cbind(new_matrix,A_indiv,B0_indiv,phase_indiv)
  colnames(new_matrix)[c(ncol(new_matrix)-2,ncol(new_matrix)-1,ncol(new_matrix))]<-c(paste0(unique_tfs[tf],".A"),paste0(unique_tfs[tf],".B0"),paste0(unique_tfs[tf],".phase"))
}
