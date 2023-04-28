library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
library(ggplot2)
library(ggpubr)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

####Sample ID----
meta<-read.csv("../../../Metadata/Meta_with_clinical_and_Summary_good_samples.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="NOR_SZBA" | sample_type=="LC_SZBU" | sample_type=="NOR_SZDW")
meta$group<-meta$Group
meta$group<-gsub("Lung cancer","Lung_cancer",meta$group)
table(meta$sample_type)

###Discovery cohort
dis_case<-subset(meta,sample_type=="LC_SZDE")
dis_control<-subset(meta,sample_type=="NOR_SZBA")
dis_ID<-c(dis_case$ID,dis_control$ID)
dis_type<-factor(c(dis_case$group,dis_control$group),levels = c("Healthy","Lung_cancer"))
ndiscovery<-length(dis_ID)
###Validation cohort
vali_1<-subset(meta,sample_type=="LC_SZBU")
vali_2<-subset(meta,sample_type=="NOR_SZDW")
vali_ID<-c(vali_1$ID,vali_2$ID)
vali_type<-factor(c(vali_1$group,vali_2$group),levels = c("Healthy","Lung_cancer"))
nvalidation<-length(vali_ID)

####goi----
goi_msRNA<-readLines("../LASSO/goi/20/msRNA name.txt")
goi_miRNA<-readLines("../Boruta/goi/20/miRNA name.txt")
goi_snRNA<-readLines("../LASSO/goi/20/snRNA name.txt")
goi_snoRNA<-readLines("../TopN/goi/snoRNA Top N selected name.txt")
goi_tsRNA<-readLines("../LASSO/goi/20/tsRNA name.txt")
goi<-c(goi_msRNA,goi_miRNA,goi_snRNA,goi_snoRNA,goi_tsRNA)
####Data for machine learning----
AllRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[goi,meta$ID]

####Single RNA ROC----
seeds<-1:100 #repeat 100 times
##msRNA###LR 
subset <- t(AllRNA[goi_msRNA,])
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
validation_x <- subset[vali_ID,]
validation_y <- vali_type
registerDoParallel(cl<-makeCluster(10))
msRNA_result<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
  Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
  prob <- predict(Ridge, validation_x, type = "response")
  pred <- prediction(prob, validation_y)
  AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  AUC_validation
}
stopCluster(cl)
msRNA_result=as.data.frame(msRNA_result)
msRNA_result$seed<-1:100
msRNA_result=msRNA_result[order(-msRNA_result[,1]),]
seed=msRNA_result[50,2]#median
set.seed(seed)
splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
training_x <- discovery_x[splitSample,]
training_y <- discovery_y[splitSample]
testing_x <- discovery_x[-splitSample,]
testing_y <- discovery_y[-splitSample]
cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
LR_prob_msRNA <- predict(Ridge, validation_x, type = "response")[,1]
LR_roc_msRNA <- plot.roc(validation_y,LR_prob_msRNA)
LR_pred_msRNA <- prediction(LR_prob_msRNA, validation_y)
LR_AUC_msRNA <- round(attr(performance(LR_pred_msRNA, "auc"), "y.values")[[1]],3)
##miRNA-LR
subset <- t(AllRNA[goi_miRNA,])
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
validation_x <- subset[vali_ID,]
validation_y <- vali_type
registerDoParallel(cl<-makeCluster(10))
miRNA_result<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
  Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
  prob <- predict(Ridge, validation_x, type = "response")
  pred <- prediction(prob, validation_y)
  AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  AUC_validation
}
stopCluster(cl)
miRNA_result=as.data.frame(miRNA_result)
miRNA_result$seed<-1:100
miRNA_result=miRNA_result[order(-miRNA_result[,1]),]
seed=miRNA_result[50,2]#median
set.seed(seed)
splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
training_x <- discovery_x[splitSample,]
training_y <- discovery_y[splitSample]
testing_x <- discovery_x[-splitSample,]
testing_y <- discovery_y[-splitSample]
cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
LR_prob_miRNA <- predict(Ridge, validation_x, type = "response")[,1]
LR_roc_miRNA <- plot.roc(validation_y,LR_prob_miRNA)
LR_pred_miRNA <- prediction(LR_prob_miRNA, validation_y)
LR_AUC_miRNA <- round(attr(performance(LR_pred_miRNA, "auc"), "y.values")[[1]],3)
##tsRNA-RF
subset <- t(AllRNA[goi_tsRNA,])
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
validation_x <- subset[vali_ID,]
validation_y <- vali_type
registerDoParallel(cl<-makeCluster(10))
tsRNA_result<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","randomForest","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  RF = randomForest(x = training_x, y = training_y)
  prob <- predict(RF, validation_x, type = "prob")[,2]
  pred <- prediction(prob, validation_y)
  AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  AUC_validation
}
stopCluster(cl)
tsRNA_result=as.data.frame(tsRNA_result)
tsRNA_result$seed<-1:100
tsRNA_result=tsRNA_result[order(-tsRNA_result[,1]),]
tsRNA_result$rank<-1:100
seed=tsRNA_result[50,2]#median
set.seed(seed)
splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
training_x <- discovery_x[splitSample,]
training_y <- discovery_y[splitSample]
testing_x <- discovery_x[-splitSample,]
testing_y <- discovery_y[-splitSample]
RF <- randomForest(x = training_x, y = training_y)
RF_prob_tsRNA <- predict(RF, validation_x, type = "prob")[,2]
RF_roc_tsRNA <- plot.roc(validation_y,RF_prob_tsRNA)
RF_pred_tsRNA <- prediction(RF_prob_tsRNA, validation_y)
RF_AUC_tsRNA <- round(attr(performance(RF_pred_tsRNA, "auc"), "y.values")[[1]],3)
##snRNA
subset <- t(AllRNA[goi_snRNA,])
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
validation_x <- subset[vali_ID,]
validation_y <- vali_type
registerDoParallel(cl<-makeCluster(10))
snRNA_result<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
  SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
  prob <- predict(SVM, validation_x, type = "prob")[,2]
  pred <- prediction(prob, validation_y)
  AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  AUC_validation
}
stopCluster(cl)
snRNA_result=as.data.frame(snRNA_result)
snRNA_result$seed<-1:100
snRNA_result=snRNA_result[order(-snRNA_result[,1]),]
seed=snRNA_result[50,2]#median
set.seed(seed)
splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
training_x <- discovery_x[splitSample,]
training_y <- discovery_y[splitSample]
testing_x <- discovery_x[-splitSample,]
testing_y <- discovery_y[-splitSample]
trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
SVM_prob_snRNA <- predict(SVM, validation_x, type = "prob")[,2]
SVM_roc_snRNA <- plot.roc(validation_y,SVM_prob_snRNA)
SVM_pred_snRNA <- prediction(SVM_prob_snRNA, validation_y)
SVM_AUC_snRNA <- round(attr(performance(SVM_pred_snRNA, "auc"), "y.values")[[1]],3)
##snoRNA
subset <- t(AllRNA[goi_snoRNA,])
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
validation_x <- subset[vali_ID,]
validation_y <- vali_type
registerDoParallel(cl<-makeCluster(10))
snoRNA_result<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
  SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
  prob <- predict(SVM, validation_x, type = "prob")[,2]
  pred <- prediction(prob, validation_y)
  AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
  AUC_validation
}
stopCluster(cl)
snoRNA_result=as.data.frame(snoRNA_result)
snoRNA_result$seed<-1:100
snoRNA_result=snoRNA_result[order(-snoRNA_result[,1]),]
snoRNA_result$rank<-1:100
seed=snoRNA_result[50,2]#median
set.seed(seed)
splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
training_x <- discovery_x[splitSample,]
training_y <- discovery_y[splitSample]
testing_x <- discovery_x[-splitSample,]
testing_y <- discovery_y[-splitSample]
trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
SVM_prob_snoRNA <- predict(SVM, validation_x, type = "prob")[,2]
SVM_roc_snoRNA <- plot.roc(validation_y,SVM_prob_snoRNA)
SVM_pred_snoRNA <- prediction(SVM_prob_snoRNA, validation_y)
SVM_AUC_snoRNA <- round(attr(performance(SVM_pred_snoRNA, "auc"), "y.values")[[1]],3)

#####PLOT----
pdf("ROC-Single cfRNA signatures in lung cancer detection.pdf", width = 6, height = 6)
plot(RF_roc_tsRNA, col="#4E62AB",lwd=2,lty=1,cex.lab=1.5, cex.axis=1.5)
plot(LR_roc_msRNA, add = TRUE, col="#FDB96A",lwd=2,lty=1)
plot(SVM_roc_snoRNA, add = TRUE, col="#87CFA4",lwd=2,lty=1)
plot(SVM_roc_snRNA, add = TRUE, col="#F0D43A",lwd=2,lty=1)
plot(LR_roc_miRNA, add = TRUE, col="#D6404E",lwd=2,lty=1)
legend(x=0.7, y=0.3, legend=paste0("miRNA   LR    AUC=0.920"),col="#D6404E", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.25, legend=paste0("snRNA   SVM AUC=",SVM_AUC_snRNA),col="#F0D43A", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.2, legend=paste0("snoRNA SVM AUC=0.881"),col="#87CFA4", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.15, legend=paste0("msRNA  LR    AUC=",LR_AUC_msRNA),col="#FDB96A", lty=1,lwd=2, bty="n",cex=1.2)
legend(x=0.7, y=0.1, legend=paste0("tsRNA    RF    AUC=",RF_AUC_tsRNA),col="#4E62AB", lty=1,lwd=2, bty="n",cex=1.2)
dev.off()
