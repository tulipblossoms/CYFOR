library(randomForest)
library(dplyr)
library(caret)

combined_list<-read.csv("C:/Users/kelly/Internship/clock_project/ML/data/final_filtered_combined_list.csv") #path to my data
full_clock_gene_names<-c("ENSPANG00000004699", "ENSPANG00000006662", "ENSPANG00000006882", "ENSPANG00000010418", "ENSPANG00000004474","ENSPANG00000026255", "ENSPANG00000013724", "ENSPANG00000006538", "ENSPANG00000008425", "ENSPANG00000023840","ENSPANG00000021704")
tissues<-unique(substr(colnames(fpkm_data[,3:ncol(fpkm_data)]),1,3))

final_correlations_oscillate<-c()
sensitivity<-c()
specificity<-c()
feature_importance<-matrix(0,ncol(x_var_matrix_full),length(unique(substr(combined_list$X,1,18))))
rownames(feature_importance)<-colnames(x_var_matrix_full)
colnames(feature_importance)<-unique(substr(combined_list$X,1,18))

ones<-sample(which(y==1))
zeros<-sample(which(y==0))
ones_folds<-split(ones, cut(seq_along(ones), 5, labels = FALSE))
zeros_folds<-split(zeros, cut(seq_along(zeros), 5, labels = FALSE))
folds <- list()
for (n in 1:5) {
  folds[[n]] <- c(ones_folds[[n]], zeros_folds[[n]])
}

#start 5-fold train
predicted_probs_all<-list()
rf_model_list<-list()
predictions_full<-rep(0,length(54464))
actuals<-rep(0,length(54464))
for(j in 1:5){
  test_inds<-sample(folds[[j]])
  train_inds<-setdiff(seq(1:nrow(new_matrix)),test_inds)
  train_x<-new_matrix[train_inds,]
  train_y<-y_b0[train_inds]
  train_data<-as.data.frame(cbind(train_y,train_x))
  test_x<-new_matrix[test_inds,]
  test_y<-y_b0[test_inds]
  colnames(train_data)[1]<-"train_y"
  rf_model<-randomForest(train_y~.,data=train_data)
  #predicted_probs<-predict(rf_model,test_x,type="prob")
  predictions<-predict(rf_model,test_x)
  predicted_probs_all[[paste0("fold",j)]]<-predictions
  rf_model_list[[paste0("fold",j)]]<-rf_model
  predictions_full[test_inds]<-as.numeric(as.character(predictions))
  actuals[test_inds]<-test_y
  conf_matrix<-confusionMatrix(as.factor(predictions),as.factor(test_y))
  sensitivity<-c(sensitivity,conf_matrix$byClass["Sensitivity"])
  specificity<-c(specificity,conf_matrix$byClass["Specificity"])

}

#draw ROC Curves & Calculate AUROC
sensitivities<-c()
specificities<-c()
for(i in seq(1,length(actual_labels),5)){
  actual_labels_roc<-rep(0,length(actual_labels))
  actual_labels_roc[1:i]<-1
  conf<-confusionMatrix(factor(actual_labels_roc,levels=c("0","1")),factor(test_y,levels=c("0","1")))
  sensitivities<-c(sensitivities,conf$byClass["Sensitivity"])
  specificities<-c(specificities,conf$byClass["Specificity"])
}
plot(1-specificities,sensitivities,type="l",main="Average ROC Curve for Model",col="red")
abline(a=0,b=1)
auc <- sum(diff(rev(1-specificities)) * (head(rev(sensitivities), -1) + tail(rev(sensitivities), -1)) / 2)
text(x = 0.6, y = 0.2, labels = paste("AUC =", round(auc, 3)), col = "black", cex = 1.2)
