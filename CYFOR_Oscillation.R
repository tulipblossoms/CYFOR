library(readr)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(patchwork)

#read in data
fpkm_data<-read.csv("C:/Users/kelly/Internship/clock_project/ML/data/GSE98965_baboon_tissue_expression_FPKM.csv.gz")
rownames(fpkm_data)<-fpkm_data$EnsemblID
gene_metadata<-cbind(fpkm_data$EnsemblID,fpkm_data$Symbol)
fpkm_data<-fpkm_data[,3:ncol(fpkm_data)]
fpkm_data_log<-log2(fpkm_data+1)
inds<-seq(1,nrow(fpkm_data_log)

#fit oscillation models
new_combined_list_nt<-rep(0,5)
combined_plot_list<-list()
for(i in inds){
table_nums<-rep(0,5)
  count<-0
  plot_list<-list()
  for(j in target_tissues){
    if(j=="IRI"){
      x<-seq(0,21,by=2)
    }else{
      x<-seq(0,23,by=2)
    }
    mini_x<-seq(0,24,by=0.01)
    y<-as.numeric(fpkm_data_log[i,grepl(j,colnames(fpkm_data_log))])
    sin_values <- sin((2 * pi / 24) * x)
    cos_values <- cos((2 * pi / 24) * x)
    fit<-lm(y~sin_values+cos_values)
    fit_constant<-lm(y~1)
    table_nums[1:3]<-fit$coefficients
    anova_test<-anova(fit_constant,fit)
    variance<-anova_test$`Sum of Sq`[2]/anova_test$Df[2]
    table_nums[4:5]<-c(anova_test$`Pr(>F)`[2],variance)
    new_combined_list_nt<-rbind(new_combined_list_nt,table_nums)
    rownames(new_combined_list_nt)[nrow(new_combined_list_nt)]<-paste0(rownames(fpkm_data_log)[i],".",gene_metadata[which(gene_metadata[,1]==rownames(fpkm_data_log)[i]),2],".",j)
    
    new_y<-fit$coefficients[1]+sin(2*pi*mini_x/24)*fit$coefficients[2]+cos(2*pi*mini_x/24)*fit$coefficients[3]
    fun<-cbind(x,y1=as.numeric((y)))
    t_hour1<-which(substr(new_combined_list_nt$X,1,18)==rownames(fpkm_data_log)[i])
    t_hour2<-which(substr(new_combined_list_nt$X,nchar(new_combined_list_nt$X)-2,nchar(new_combined_list_nt$X))==j)
    t_ind<-intersect(t_hour1,t_hour2)
    t_hour<-new_combined_list_nt$hour[intersect(t_hour1,t_hour2)]
    color<-ifelse(new_combined_list_nt$p.val[t_ind]<0.2,"darkred","darkblue")
    p_val<-new_combined_list_nt$p.val[t_ind]
    circadian<-ifelse(new_combined_list_nt$p.val[t_ind]<0.2,"C","NC")
    
    p <- ggplot() +geom_line(data = as.data.frame(cbind(mini_x, new_y)), aes(x = mini_x, y = new_y), color = "pink", size = 2)+labs(title=circadian,x = "Zeitgeber Time (ZT) ",y = "Expression")+theme(plot.title = element_text(size = 3, hjust = 0.5))+scale_x_continuous(breaks = seq(0, 23, by = 2), labels = seq(0, 23, by = 2))+geom_point(data = as.data.frame(fun), aes(x = x, y = y1),color="red", size = 3.5)+scale_color_manual(name = "Annotations",values = c("Raw Gene Exp" = "red", "Peak Time" = "black"))+theme_minimal()+  theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),plot.title=element_text(size=30,color=color,hjust=0.5),
      axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15)) 
    
    if(p_val<0.2){
      p<-p+geom_vline(xintercept = t_hour,color = "black", linetype = "dashed",linewidth=1.5)
    }
                                                                                                                                                                                                                         
    count<-count+1
    plot_list[[count]]<-p
    
  }

  combined_plot<- Reduce(`+`, plot_list) + plot_layout(ncol = 3)
  combined_plot_list[[i]]<-combined_plot
}

#finding phase
colnames(new_combined_list_nt)<-c("X","B0","B1","B2","p.val","variance")
new_combined_list_nt<-new_combined_list_nt[2:nrow(new_combined_list_nt),]
new_combined_list_nt<-as.data.frame(new_combined_list_nt)
new_combined_list_nt$A_val<-sqrt(new_combined_list_nt$B1^2+new_combined_list_nt$B2^2)

new_combined_list_nt$new_cos<-acos(new_combined_list_nt$B2/new_combined_list_nt$A_val)
new_combined_list_nt$new_sin<-asin(new_combined_list_nt$B1/(new_combined_list_nt$A_val))

calculate_angles <- function(acos_value, asin_value) {
  if(is.na(acos_value)||is.na(asin_value)){
    angle<-NaN
  }else if (asin_value < 0){
    angle <- 2 * pi + asin_value  
  } else {
    angle <- acos_value 
  }
  return(angle)
}
new_combined_list_nt$new_final_phase <- mapply(calculate_angles, new_combined_list_nt$new_cos, new_combined_list_nt$new_sin)
#new_combined_list$new_final_phase<-new_combined_list$new_final_phase-(2*pi)
new_combined_list_nt$hour_sincos<-(new_combined_list_nt$new_final_phase*24)/(2*pi)

#TAN REINFORCEMENTS
new_combined_list_nt$new_final_tan <- atan2(new_combined_list_nt$B1, new_combined_list_nt$B2)
new_combined_list_nt$new_final_tan <- ifelse(new_combined_list_nt$new_final_tan < 0, 
                                            2 * pi + new_combined_list_nt$new_final_tan, 
                                            new_combined_list_nt$new_final_tan)
new_combined_list_nt$hour <- (new_combined_list_nt$new_final_tan * 24) / (2 * pi)

new_combined_list_nt$p.val<-p.adjust(new_combined_list_nt$p.val,method="fdr")
fpkm_significant<-new_combined_list_nt[which(new_combined_list_nt$p.val<0.2),]
