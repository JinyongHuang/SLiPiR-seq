library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#library(ComplexHeatmap)
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

####Data for machine learning----
AllRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[,meta$ID]

####goi----
goi_mRNA<-readLines("../LASSO_10/goi/mRNA name.txt")
goi_miRNA<-readLines("../LASSO_10/goi/miRNA name.txt")
goi_snRNA<-readLines("../LASSO_10/goi/snRNA name.txt")
goi_snoRNA<-readLines("../LASSO_10/goi/snoRNA name.txt")
goi_tsRNA<-readLines("../LASSO_10/goi/tsRNA name.txt")
# goi_rsRNA<-readLines("../LASSO_10/goi/rsRNA name.txt")
# goi_ysRNA<-readLines("../LASSO_10/goi/ysRNA name.txt")

####Different combination----
goi_combn=list(goi_mRNA=goi_mRNA,goi_miRNA=goi_miRNA,goi_snRNA=goi_snRNA,goi_snoRNA=goi_snoRNA,goi_tsRNA=goi_tsRNA)
combn1<-data.frame(t(combn(names(goi_combn),1)))
combn1$name<-substr(combn1[,1],5,nchar(combn1[,1])-3)
combn2<-data.frame(t(combn(names(goi_combn),2)))
combn2$name<-paste0(substr(combn2[,1],5,nchar(combn2[,1])-3),"+",substr(combn2[,2],5,nchar(combn2[,2])-3))
combn3<-data.frame(t(combn(names(goi_combn),3)))
combn3$name<-paste0(substr(combn3[,1],5,nchar(combn3[,1])-3),"+",substr(combn3[,2],5,nchar(combn3[,2])-3),"+",substr(combn3[,3],5,nchar(combn3[,3])-3))
combn4<-data.frame(t(combn(names(goi_combn),4)))
combn4$name<-paste0(substr(combn4[,1],5,nchar(combn4[,1])-3),"+",substr(combn4[,2],5,nchar(combn4[,2])-3),"+",
                    substr(combn4[,3],5,nchar(combn4[,3])-3),"+",substr(combn4[,4],5,nchar(combn4[,4])-3))
combn5<-data.frame(t(combn(names(goi_combn),5)))
combn5$name<-paste0(substr(combn5[,1],5,nchar(combn5[,1])-3),"+",substr(combn5[,2],5,nchar(combn5[,2])-3),"+",
                    substr(combn5[,3],5,nchar(combn5[,3])-3),"+",substr(combn5[,4],5,nchar(combn5[,4])-3),"+",substr(combn5[,5],5,nchar(combn5[,5])-3))
# combn6<-data.frame(t(combn(names(goi_combn),6)))
# combn6$name<-paste0(substr(combn6[,1],5,nchar(combn6[,1])-3),"+",substr(combn6[,2],5,nchar(combn6[,2])-3),"+",substr(combn6[,3],5,nchar(combn6[,3])-3),"+",
#                     substr(combn6[,4],5,nchar(combn6[,4])-3),"+",substr(combn6[,5],5,nchar(combn6[,5])-3),"+",substr(combn6[,6],5,nchar(combn6[,6])-3))
# combn7<-data.frame(t(combn(names(goi_combn),7)))
# combn7$name<-paste0(substr(combn7[,1],5,nchar(combn7[,1])-3),"+",substr(combn7[,2],5,nchar(combn7[,2])-3),"+",substr(combn7[,3],5,nchar(combn7[,3])-3),"+",
#                     substr(combn7[,4],5,nchar(combn7[,4])-3),"+",substr(combn7[,5],5,nchar(combn7[,5])-3),"+",substr(combn7[,6],5,nchar(combn7[,6])-3),"+",substr(combn7[,7],5,nchar(combn7[,6])-3))
Combine<-list()
for (i in 1:nrow(combn1)) {
  Combine[[combn1$name[i]]]=c(goi_combn[[combn1[,1][i]]])
}
for (i in 1:nrow(combn2)) {
  Combine[[combn2$name[i]]]=c(goi_combn[[combn2[,1][i]]],goi_combn[[combn2[,2][i]]])
}
for (i in 1:nrow(combn3)) {
  Combine[[combn3$name[i]]]=c(goi_combn[[combn3[,1][i]]],goi_combn[[combn3[,2][i]]],goi_combn[[combn3[,3][i]]])
}
for (i in 1:nrow(combn4)) {
  Combine[[combn4$name[i]]]=c(goi_combn[[combn4[,1][i]]],goi_combn[[combn4[,2][i]]],goi_combn[[combn4[,3][i]]],goi_combn[[combn4[,4][i]]])
}
for (i in 1:nrow(combn5)) {
  Combine[[combn5$name[i]]]=c(goi_combn[[combn5[,1][i]]],goi_combn[[combn5[,2][i]]],goi_combn[[combn5[,3][i]]],goi_combn[[combn5[,4][i]]],goi_combn[[combn5[,5][i]]])
}
# for (i in 1:nrow(combn6)) {
#   Combine[[combn6$name[i]]]=c(goi_combn[[combn6[,1][i]]],goi_combn[[combn6[,2][i]]],goi_combn[[combn6[,3][i]]],goi_combn[[combn6[,4][i]]],goi_combn[[combn6[,5][i]]],goi_combn[[combn6[,6][i]]])
# }
# for (i in 1:nrow(combn7)) {
#   Combine[[combn7$name[i]]]=c(goi_combn[[combn7[,1][i]]],goi_combn[[combn7[,2][i]]],goi_combn[[combn7[,3][i]]],goi_combn[[combn7[,4][i]]],goi_combn[[combn7[,5][i]]],goi_combn[[combn7[,6][i]]],goi_combn[[combn7[,7][i]]])
# }
####Machine learning----
seeds<-1:100 #repeat 100 times
Combination<-list()
for(i in names(Combine)){
  print(paste("Start working on",i))
  top_rpm<-AllRNA[Combine[[i]],]
  top_rpm <- t(top_rpm)
  discovery_x <- top_rpm[dis_ID,]
  discovery_y <- dis_type
  validation_x <- top_rpm[vali_ID,]
  validation_y <- vali_type
  ###LR logistic regression
  registerDoParallel(cl<-makeCluster(10))
  result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
    set.seed(seed) #for loop
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
    LR <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
    ##Training set
    pred.class <- predict(LR, training_x, type = "class")
    ACC_training<-round(100* mean(pred.class[,1] == training_y),2)
    prob <- predict(LR, training_x, type = "response")
    pred <- prediction(prob, training_y)
    AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(training_y,as.numeric(prob))
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_training<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    ##Held-out testing set
    pred.class <- predict(LR, testing_x, type = "class")
    ACC_testing<-round(100* mean(pred.class[,1] == testing_y),2)
    prob <- predict(LR, testing_x, type = "response")
    pred <- prediction(prob,testing_y)
    AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(testing_y,as.numeric(prob))
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_testing<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    ##Validation set
    pred.class <- predict(LR, validation_x, type = "class")
    ACC_validation<-round(100* mean(pred.class[,1] == validation_y),2)
    prob <- predict(LR, validation_x, type = "response")
    pred <- prediction(prob, validation_y)
    AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(validation_y,as.numeric(prob))
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_validation<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    c(ACC_training,ACC_testing,ACC_validation,AUC_training,AUC_testing,AUC_validation,SenAt1Spe_training,SenAt1Spe_testing,SenAt1Spe_validation)
  }
  stopCluster(cl)
  colnames(result_seed)= c('ACC_training','ACC_testing','ACC_validation','AUC_training','AUC_testing','AUC_validation','SenAt1Spe_training','SenAt1Spe_testing','SenAt1Spe_validation')
  LR_df<-result_seed
  
  ###Random forest
  registerDoParallel(cl<-makeCluster(10))
  result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","randomForest","ROCR","pROC")) %dopar% {
    set.seed(seed) #for loop
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    RF = randomForest(x = training_x, y = training_y)
    ##Training set
    pred.class <- predict(RF, training_x, type = "class")
    ACC_training<-round(100* mean(pred.class == training_y),2)
    prob <- predict(RF, training_x, type = "prob")
    pred <- prediction(prob[,2], training_y)
    AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(training_y,prob[,2])
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_training<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    ##Held-out testing set
    pred.class <- predict(RF, testing_x, type = "class")
    ACC_testing<-round(100* mean(pred.class == testing_y),2)
    prob <- predict(RF, testing_x, type = "prob")
    pred <- prediction(prob[,2], testing_y)
    AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(testing_y,prob[,2])
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_testing<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    ##Validation set
    pred.class <- predict(RF, validation_x, type = "class")
    ACC_validation<-round(100* mean(pred.class == validation_y),2)
    prob <- predict(RF, validation_x, type = "prob")
    pred <- prediction(prob[,2], validation_y)
    AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(validation_y,prob[,2])
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_validation<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    c(ACC_training,ACC_testing,ACC_validation,AUC_training,AUC_testing,AUC_validation,SenAt1Spe_training,SenAt1Spe_testing,SenAt1Spe_validation)
  }
  stopCluster(cl)
  colnames(result_seed)= c('ACC_training','ACC_testing','ACC_validation','AUC_training','AUC_testing','AUC_validation','SenAt1Spe_training','SenAt1Spe_testing','SenAt1Spe_validation')
  RF_df<-result_seed
  
  ###Support vector machine
  registerDoParallel(cl<-makeCluster(10))
  result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","ROCR","pROC")) %dopar% {
    set.seed(seed) #for loop
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
    SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
    ##Training set
    pred.class <- predict(SVM, training_x, type = "raw")
    ACC_training<-round(100* mean(pred.class == training_y),2)
    prob <- predict(SVM, training_x, type = "prob")
    pred <- prediction(prob[,2], training_y)
    AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(training_y,prob[,2])
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_training<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    ##Held-out testing set
    pred.class <- predict(SVM, testing_x, type = "raw")
    ACC_testing<-round(100* mean(pred.class == testing_y),2)
    prob <- predict(SVM, testing_x, type = "prob")
    pred <- prediction(prob[,2], testing_y)
    AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(testing_y,prob[,2])
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_testing<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    ##Validation set
    pred.class <- predict(SVM, validation_x, type = "raw")
    ACC_validation<-round(100* mean(pred.class == validation_y),2)
    prob <- predict(SVM, validation_x, type = "prob")
    pred <- prediction(prob[,2], validation_y)
    AUC_validation <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
    roc<-roc(validation_y,prob[,2])
    coord_list <- coords(roc, x = "all")
    SenAt1Spe<-subset(coord_list,specificity==1)#sensitivity when specificity=1
    SenAt1Spe_validation<-SenAt1Spe[order(-SenAt1Spe$sensitivity),][1,3]
    c(ACC_training,ACC_testing,ACC_validation,AUC_training,AUC_testing,AUC_validation,SenAt1Spe_training,SenAt1Spe_testing,SenAt1Spe_validation)
  }
  stopCluster(cl)
  colnames(result_seed)= c('ACC_training','ACC_testing','ACC_validation','AUC_training','AUC_testing','AUC_validation','SenAt1Spe_training','SenAt1Spe_testing','SenAt1Spe_validation')
  SVM_df<-result_seed
  
  Combination[[i]]<-data.frame(rbind(LR_df,RF_df,SVM_df))
}

dir.create(paste0(getwd(),"/AUC plot all in one"), showWarnings = FALSE)
for (i in names(Combine)) {  
  n=length(Combine[[i]])
  df<-Combination[[i]][,5:6]
  df$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  df<-melt(df)
  df$variable<-factor(c(rep("Testing",300),rep("Validation",300)),levels = c("Training","Testing","Validation"))
  ggplot(df, aes(x=variable, y= value,colour=variable))+geom_boxplot()+  geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(~method, scale="free")+
    labs(x="", y = "AUC in 100 iterations")+ theme_classic()+scale_colour_brewer(palette="Set1")+
    theme(axis.text.y = element_text(colour = "black",size = 14), axis.text.x = element_text(colour = "black",size = 14),
          strip.text.x = element_text(colour = "black",size = 14,face = "bold"), legend.position = "none", 
          axis.title = element_text(size = 16), plot.title = element_text(size = 16,hjust = 0.5,face = "bold"), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  ggsave(paste0("AUC plot all in one/Combination of ", i, " N=",n," - AUC.pdf"), width = 7, height = 7)
  ggsave(paste0("AUC plot all in one/Combination of ", i, " N=",n," - AUC.png"), width = 7, height = 7, bg="white")
}

### boxplot of all combinations----
#AUC_testing
All_AUC<-lapply(Combination, as.data.frame)
All_AUC<-as.data.frame(lapply(All_AUC, function(x) x[,"AUC_testing"]))
colnames(All_AUC)<-gsub("\\.","+",colnames(All_AUC))
All_AUC_median<-apply(All_AUC, 2, median)
All_AUC_median<-sort(All_AUC_median,decreasing = T)
All_AUC$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
df<-melt(All_AUC)
df$variable<-factor(df$variable,levels = names(All_AUC_median))
ggplot(df, aes(x=variable, y= value))+geom_boxplot(outlier.shape = NA)+  geom_jitter(aes(color=method),shape=16, position=position_jitter(0.2),alpha=0.8,size=0.8) +
  labs(x="", y = "AUC in the test set")+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+
  theme(axis.text.y = element_text(colour = "black",size = 12), axis.text.x = element_text(colour = "black",size = 12,angle = 60,hjust=1),
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = c(0.03,0.15), strip.text = element_blank(),
        axis.title = element_text(size = 12), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(paste0("AUC plot AUC_testing all combinations.pdf"), width = 15, height = 5)
ggsave(paste0("AUC plot AUC_testing all combinations.png"), width = 15, height = 5, bg="white")
#AUC_validation
All_AUC<-lapply(Combination_sub, as.data.frame)
All_AUC<-as.data.frame(lapply(All_AUC, function(x) x[,"AUC_validation"]))
colnames(All_AUC)<-gsub("\\.","+",colnames(All_AUC))
All_AUC_median<-apply(All_AUC, 2, median)
All_AUC_median<-sort(All_AUC_median,decreasing = T)
All_AUC$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
df<-melt(All_AUC)
df$variable<-factor(df$variable,levels = names(All_AUC_median))
ggplot(df, aes(x=variable, y= value))+geom_boxplot(outlier.shape = NA)+  geom_jitter(aes(color=method),shape=16, position=position_jitter(0.2),alpha=0.8,size=0.8) +
  labs(x="", y = "AUC in the validation cohort")+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+
  theme(axis.text.y = element_text(colour = "black",size = 12), axis.text.x = element_text(colour = "black",size = 12,angle = 60,hjust=1),
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = c(0.03,0.15),strip.text = element_blank(),
        axis.title = element_text(size = 12), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(paste0("AUC plot AUC_validation all combinations.pdf"), width = 15, height = 5)
ggsave(paste0("AUC plot AUC_validation all combinations.png"), width = 15, height = 5, bg="white")
#AUC_testing--LR 
#AUC_validation--LR 
All_AUC_test<-lapply(Combination_sub, as.data.frame)
All_AUC_test<-as.data.frame(lapply(All_AUC_test, function(x) x[,"AUC_testing"]))
colnames(All_AUC_test)<-gsub("\\.","+",colnames(All_AUC_test))
All_AUC_test<-All_AUC_test[1:100,]#LR
All_AUC_vali<-lapply(Combination_sub, as.data.frame)
All_AUC_vali<-as.data.frame(lapply(All_AUC_vali, function(x) x[,"AUC_validation"]))
colnames(All_AUC_vali)<-gsub("\\.","+",colnames(All_AUC_vali))
All_AUC_vali<-All_AUC_vali[1:100,]#LR
All_AUC<-rbind(All_AUC_test,All_AUC_vali)
All_AUC_median<-apply(All_AUC, 2, median)
All_AUC_median<-sort(All_AUC_median,decreasing = T)
All_AUC$cohort<-c(rep("Discovery-test",100),rep("Validation",100))
df<-melt(All_AUC)
df$variable<-factor(df$variable,levels = names(All_AUC_median))
df<-df[df$variable %in% names(All_AUC_median),]
ggplot(df, aes(x=variable, y= value))+geom_boxplot(outlier.shape = NA)+  geom_jitter(aes(color=cohort),shape=16, position=position_jitter(0.2),size=0.8) +
  labs(x="", y = "AUC in LR models")+ theme_classic()+scale_colour_manual(values=c("#4E62AB","#D6404E"))+
  theme(axis.text.y = element_text(colour = "black",size = 12), axis.text.x = element_text(colour = "black",size = 12,angle = 60,hjust=1),
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = c(0.1,0.15),strip.text = element_blank(), legend.background = element_blank(),
        axis.title = element_text(size = 12), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(paste0("AUC plot LR AUC all combinations.pdf"), width = 8, height = 5)
ggsave(paste0("AUC plot LR AUC all combinations.png"), width = 8, height = 5, bg="white")


###summary----
for (i in names(Combine)) {
  df<-Combination[[i]][,4:6]
  df$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  summary<-aggregate(df[,1:3], by=list(method=df$method), FUN = median)
  summary$Name=i
  if (i==names(Combine)[1]){AUCSummary=summary}
  if (i!=names(Combine)[1]){AUCSummary=rbind(AUCSummary,summary)}
}
write.csv(AUCSummary,"Combination AUC summary.csv")

####Risk score boxplot----
dir.create(paste0(getwd(),"/Risk score boxplot"), showWarnings = FALSE)
seeds<-1:100 #repeat 100 times
####LR Risk score
dir.create(paste0(getwd(),"/Risk score boxplot/LR"), showWarnings = FALSE)
prob_discovery_LR_list<-list()
prob_validation_LR_list<-list()
for(i in names(Combine)){
  print(paste("Start working on",i))
  top_rpm<-AllRNA[Combine[[i]],]
  top_rpm <- t(top_rpm)
  discovery_x <- top_rpm[dis_ID,]
  discovery_y <- dis_type
  validation_x <- top_rpm[vali_ID,]
  validation_y <- vali_type
  registerDoParallel(cl<-makeCluster(10))
  prob_result<-foreach(seed=seeds, .combine="cbind", .packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
    set.seed(seed)
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
    Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
    ##Discovery cohort
    prob_discovery <- predict(Ridge, discovery_x, type = "response")
    ##Validation cohort
    prob_validation <- predict(Ridge, validation_x, type = "response")
    c(prob_discovery,prob_validation)
  }
  stopCluster(cl)
  prob_discovery_LR<-prob_result[1:length(discovery_y),]
  prob_discovery_LR_list[[i]]<-data.frame(prob_discovery_LR)
  prob_validation_LR<-prob_result[(length(discovery_y)+1):(length(discovery_y)+length(validation_y)),]
  prob_validation_LR_list[[i]]<-data.frame(prob_validation_LR)
}
####boxplot 
for(i in names(Combine)){
  prob_discovery_LR<-prob_discovery_LR_list[[i]]
  prob_validation_LR<-prob_validation_LR_list[[i]]
  
  df_discovery<-data.frame(id=rownames(discovery_x),prob_discovery_LR,row.names = rownames(discovery_x))
  df_discovery_median<-apply(df_discovery[,-1], 1, median)
  df_discovery_median<-sort(df_discovery_median)
  df_discovery_plot<-melt(df_discovery,"id")
  df_discovery_plot$id<-factor(df_discovery_plot$id,levels = names(df_discovery_median))
  df_discovery_plot$actual=substr(df_discovery_plot$id,1,2)
  df_discovery_plot$actual=gsub("LC","Lung cancer",df_discovery_plot$actual)
  df_discovery_plot$actual=gsub("NO","Cancer-free",df_discovery_plot$actual)
  df_discovery_plot$actual=factor(df_discovery_plot$actual,levels = c("Lung cancer","Cancer-free"))
  discovery_plot<-ggplot(df_discovery_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.2,width=0.85)+ scale_colour_manual(values=c("#D6404E","#4E62AB"))+
    labs(x=paste0("Study subjects in the discovery cohort (N=",length(discovery_y),")"), y="Lung cancer risk scores")+  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
    theme_classic()+theme(axis.text.x = element_blank(), axis.text.y = element_text(colour = "black"), text = element_text(size = 8), 
                          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.15),
                          legend.position = c(0.10, 0.85), legend.title = element_blank(),legend.background =element_blank() )
  discovery_plot
  ggsave(paste0("Risk score boxplot/LR/","Combination of ", i,"-discovery-logistic regression.pdf"), width = 7, height = 2.5)
  ggsave(paste0("Risk score boxplot/LR/","Combination of ", i,"-discovery-logistic regression.png"), width = 7, height = 2.5, bg="white")
  
  df_validation<-data.frame(id=rownames(validation_x),prob_validation_LR,row.names = rownames(validation_x))
  df_validation_median<-apply(df_validation[,-1], 1, median)
  df_validation_median<-sort(df_validation_median)
  df_validation_plot<-melt(df_validation,"id")
  df_validation_plot$id<-factor(df_validation_plot$id,levels = names(df_validation_median))
  df_validation_plot$actual=substr(df_validation_plot$id,1,2)
  df_validation_plot$actual=gsub("LC","Lung cancer",df_validation_plot$actual)
  df_validation_plot$actual=gsub("NO","Cancer-free",df_validation_plot$actual)
  df_validation_plot$actual=factor(df_validation_plot$actual,levels = c("Lung cancer","Cancer-free"))
  validation_plot<-ggplot(df_validation_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.2,width=0.85)+ scale_colour_manual(values=c("#D6404E","#4E62AB"))+
    labs(x=paste0("Study subjects in the validation cohort (N=",length(validation_y),")"), y="Lung cancer risk score")+  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
    theme_classic()+theme(axis.text.x = element_blank(), axis.text.y = element_text(colour = "black"), text = element_text(size = 14), 
                          panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
                          legend.position = c(0.2, 0.87), legend.title = element_blank())
  validation_plot
  ggsave(paste0("Risk score boxplot/LR/","Combination of ", i,"-validation-logistic regression.pdf"), width = 4.8, height = 4.2)
  ggsave(paste0("Risk score boxplot/LR/","Combination of ", i,"-validation-logistic regression.png"), width = 4.8, height = 4.2, bg="white")
}  
####RF Risk score
dir.create(paste0(getwd(),"/Risk score boxplot/RF"), showWarnings = FALSE)
prob_discovery_RF_list<-list()
prob_validation_RF_list<-list()
for(i in names(Combine)){
  print(paste("Start working on",i))
  top_rpm<-AllRNA[Combine[[i]],]
  top_rpm <- t(top_rpm)
  discovery_x <- top_rpm[dis_ID,]
  discovery_y <- dis_type
  validation_x <- top_rpm[vali_ID,]
  validation_y <- vali_type  
  registerDoParallel(cl<-makeCluster(10))
  prob_result<-foreach(seed=seeds, .combine="cbind", .packages=c("tidyverse","caret","randomForest","ROCR","pROC")) %dopar% {
    set.seed(seed)
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    RF = randomForest(x = training_x, y = training_y)
    ##Discovery cohort
    prob_discovery <- predict(RF, discovery_x, type = "prob")[,2]
    ##Validation cohort
    prob_validation <- predict(RF, validation_x, type = "prob")[,2]
    c(prob_discovery,prob_validation)
  }
  stopCluster(cl)
  prob_discovery_RF<-prob_result[1:length(discovery_y),]
  prob_discovery_RF_list[[i]]<-data.frame(prob_discovery_RF)
  prob_validation_RF<-prob_result[(length(discovery_y)+1):(length(discovery_y)+length(validation_y)),]
  prob_validation_RF_list[[i]]<-data.frame(prob_validation_RF)
}
  ####boxplot
for(i in names(Combine)){
  prob_discovery_RF<-prob_discovery_RF_list[[i]]
  prob_validation_RF<-prob_validation_RF_list[[i]]
#   df_discovery<-data.frame(id=rownames(discovery_x),prob_discovery_RF,row.names = rownames(discovery_x))
#   df_discovery_median<-apply(df_discovery[,-1], 1, median)
#   df_discovery_median<-sort(df_discovery_median)
#   df_discovery_plot<-melt(df_discovery,"id")
#   df_discovery_plot$id<-factor(df_discovery_plot$id,levels = names(df_discovery_median))
#   df_discovery_plot$actual=substr(df_discovery_plot$id,1,3)
#   df_discovery_plot$actual=gsub("_","",df_discovery_plot$actual)
#   df_discovery_plot$actual=gsub("LC","Lung cancer",df_discovery_plot$actual)
#   df_discovery_plot$actual=gsub("NC","Cancer-free",df_discovery_plot$actual)
#   df_discovery_plot$actual=factor(df_discovery_plot$actual,levels = c("Lung cancer","Cancer-free"))
#   discovery_plot<-ggplot(df_discovery_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.6)+ scale_color_brewer(type = "Qualitative", palette = "Set1")+
#     labs(x=paste0("Study subjects in discovery cohort (N=",length(discovery_y),")"), y="Lung cancer risk score", title=paste0("Combination of ", i,"-discovery-random forest"))+ geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
#     theme_classic()+theme(axis.text.x = element_blank(), text = element_text(size = 14), panel.border = element_rect(colour = "black", fill=NA, linewidth=1), legend.position = c(0.08, 0.77), legend.title = element_blank())

#   ggsave(paste0("Risk score boxplot/RF/","Combination of ", i,"-discovery-random forest.pdf"), width = 18, height = 6, scale=0.8)
#   ggsave(paste0("Risk score boxplot/RF/","Combination of ", i,"-discovery-random forest.png"), width = 18, height = 6, bg="white", scale=0.8)

  df_validation<-data.frame(id=rownames(validation_x),prob_validation_RF,row.names = rownames(validation_x))
  df_validation_median<-apply(df_validation[,-1], 1, median)
  df_validation_median<-sort(df_validation_median)
  df_validation_plot<-melt(df_validation,"id")
  df_validation_plot$id<-factor(df_validation_plot$id,levels = names(df_validation_median))
  df_validation_plot$actual=substr(df_validation_plot$id,1,2)
  df_validation_plot$actual=gsub("LC","Lung cancer",df_validation_plot$actual)
  df_validation_plot$actual=gsub("NO","Cancer-free",df_validation_plot$actual)
  df_validation_plot$actual=factor(df_validation_plot$actual,levels = c("Lung cancer","Cancer-free"))
  validation_plot<-ggplot(df_validation_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.2,width=0.85)+ scale_colour_manual(values=c("#D6404E","#4E62AB"))+
    labs(x=paste0("Study subjects in validation cohort (N=",length(validation_y),")"), y="Lung cancer risk score")+  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
    theme_classic()+theme(axis.text.x = element_blank(), text = element_text(size = 14), panel.border = element_rect(colour = "black", fill=NA, linewidth=1), legend.position = c(0.2, 0.87), legend.title = element_blank())
  validation_plot
  ggsave(paste0("Risk score boxplot/RF/","Combination of ", i,"-validation-random forest.pdf"), width = 4.8, height = 4.2)
  ggsave(paste0("Risk score boxplot/RF/","Combination of ", i,"-validation-random forest.png"), width = 4.8, height = 4.2, bg="white")
  }
  ####SVM Risk score
dir.create(paste0(getwd(),"/Risk score boxplot/SVM"), showWarnings = FALSE)
prob_discovery_SVM_list<-list()
prob_validation_SVM_list<-list()
for(i in names(Combine)){
  print(paste("Start working on",i))
  top_rpm<-AllRNA[Combine[[i]],]
  top_rpm <- t(top_rpm)
  discovery_x <- top_rpm[dis_ID,]
  discovery_y <- dis_type
  validation_x <- top_rpm[vali_ID,]
  validation_y <- vali_type  
  registerDoParallel(cl<-makeCluster(10))
  prob_result<-foreach(seed=seeds, .combine="cbind", .packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
    set.seed(seed)
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
    SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
    ##Discovery cohort
    prob_discovery <- predict(SVM, discovery_x, type = "prob")[,2]
    ##Validation cohort
    prob_validation <- predict(SVM, validation_x, type = "prob")[,2]
    c(prob_discovery,prob_validation)
  }
  stopCluster(cl)
  prob_discovery_SVM<-prob_result[1:length(discovery_y),]
  prob_discovery_SVM_list[[i]]<-data.frame(prob_discovery_SVM)
  prob_validation_SVM<-prob_result[(length(discovery_y)+1):(length(discovery_y)+length(validation_y)),]
  prob_validation_SVM_list[[i]]<-data.frame(prob_validation_SVM)
}
  ####boxplot
for(i in names(Combine)){
  prob_discovery_SVM<-prob_discovery_SVM_list[[i]]
  prob_validation_SVM<-prob_validation_SVM_list[[i]]
#   df_discovery<-data.frame(id=rownames(discovery_x),prob_discovery_SVM,row.names = rownames(discovery_x))
#   df_discovery_median<-apply(df_discovery[,-1], 1, median)
#   df_discovery_median<-sort(df_discovery_median)
#   df_discovery_plot<-melt(df_discovery,"id")
#   df_discovery_plot$id<-factor(df_discovery_plot$id,levels = names(df_discovery_median))
#   df_discovery_plot$actual=substr(df_discovery_plot$id,1,3)
#   df_discovery_plot$actual=gsub("_","",df_discovery_plot$actual)
#   df_discovery_plot$actual=gsub("LC","Lung cancer",df_discovery_plot$actual)
#   df_discovery_plot$actual=gsub("NC","Cancer-free",df_discovery_plot$actual)
#   df_discovery_plot$actual=factor(df_discovery_plot$actual,levels = c("Lung cancer","Cancer-free"))
#   discovery_plot<-ggplot(df_discovery_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.6)+ scale_color_brewer(type = "Qualitative", palette = "Set1")+
#     labs(x=paste0("Study subjects in discovery cohort (N=",length(discovery_y),")"), y="Lung cancer risk score", title=paste0("Combination of ", i,"-discovery-support vector machine"))+ geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
#     theme_classic()+theme(axis.text.x = element_blank(), text = element_text(size = 14), panel.border = element_rect(colour = "black", fill=NA, linewidth=1), legend.position = c(0.08, 0.77), legend.title = element_blank())
#   discovery_plot
#   ggsave(paste0("Risk score boxplot/SVM/","Combination of ", i,"-discovery-support vector machine.pdf"), width = 18, height = 6, scale=0.8)
#   ggsave(paste0("Risk score boxplot/SVM/","Combination of ", i,"-discovery-support vector machine.png"), width = 18, height = 6, bg="white", scale=0.8)

  df_validation<-data.frame(id=rownames(validation_x),prob_validation_SVM,row.names = rownames(validation_x))
  df_validation_median<-apply(df_validation[,-1], 1, median)
  df_validation_median<-sort(df_validation_median)
  df_validation_plot<-melt(df_validation,"id")
  df_validation_plot$id<-factor(df_validation_plot$id,levels = names(df_validation_median))
  df_validation_plot$actual=substr(df_validation_plot$id,1,2)
  df_validation_plot$actual=gsub("LC","Lung cancer",df_validation_plot$actual)
  df_validation_plot$actual=gsub("NO","Cancer-free",df_validation_plot$actual)
  df_validation_plot$actual=factor(df_validation_plot$actual,levels = c("Lung cancer","Cancer-free"))
  validation_plot<-ggplot(df_validation_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.2,width=0.85)+ scale_colour_manual(values=c("#D6404E","#4E62AB"))+
    labs(x=paste0("Study subjects in validation cohort (N=",length(validation_y),")"), y="Lung cancer risk score")+ geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
    theme_classic()+theme(axis.text.x = element_blank(), text = element_text(size = 14), panel.border = element_rect(colour = "black", fill=NA, linewidth=1), legend.position = c(0.2, 0.87), legend.title = element_blank())
  validation_plot
  ggsave(paste0("Risk score boxplot/SVM/","Combination of ", i,"-validation-support vector machine.pdf"), width = 4.8, height = 4.2)
  ggsave(paste0("Risk score boxplot/SVM/","Combination of ", i,"-validation-support vector machine.png"), width = 4.8, height = 4.2, bg="white")
}

####Sensitivity and specificity----
result<-matrix(ncol=12,nrow = length(Combine))
colnames(result)<-c("LR_dis_sensitivity","LR_dis_specificity","LR_vali_sensitivity","LR_vali_specificity",
                    "RF_dis_sensitivity","RF_dis_specificity","RF_vali_sensitivity","RF_vali_specificity",
                    "SVM_dis_sensitivity","SVM_dis_specificity","SVM_vali_sensitivity","SVM_vali_specificity")
rownames(result)<-names(Combine)
for(i in names(Combine)){
  temp<-prob_discovery_LR_list[[i]]
  temp<-data.frame(apply(temp, 1, median))
  rownames(temp)<-rownames(discovery_x)
  temp$Actual<-NA
  temp$Actual[grep("LC",rownames(temp))]<-"LC"
  temp$Actual[grep("NOR",rownames(temp))]<-"NOR"
  temp$Predicted<-ifelse(temp$apply.temp..1..median.>=0.5,"LC","NOR")
  temp1<-temp[temp$Predicted=="LC",]
  result[i,"LR_dis_specificity"]<-mean(temp1$Predicted==temp1$Actual)
  temp2<-length(which(temp1$Actual=="LC"))
  result[i,"LR_dis_sensitivity"]<-temp2/139
  
  temp<-prob_validation_LR_list[[i]]
  temp<-data.frame(apply(temp, 1, median))
  rownames(temp)<-rownames(validation_x)
  temp$Actual<-NA
  temp$Actual[grep("LC",rownames(temp))]<-"LC"
  temp$Actual[grep("NOR",rownames(temp))]<-"NOR"
  temp$Predicted<-ifelse(temp$apply.temp..1..median.>=0.5,"LC","NOR")
  temp1<-temp[temp$Predicted=="LC",]
  result[i,"LR_vali_specificity"]<-mean(temp1$Predicted==temp1$Actual)
  temp2<-length(which(temp1$Actual=="LC"))
  result[i,"LR_vali_sensitivity"]<-temp2/26
  
  temp<-prob_discovery_RF_list[[i]]
  temp<-data.frame(apply(temp, 1, median))
  rownames(temp)<-rownames(discovery_x)
  temp$Actual<-NA
  temp$Actual[grep("LC",rownames(temp))]<-"LC"
  temp$Actual[grep("NOR",rownames(temp))]<-"NOR"
  temp$Predicted<-ifelse(temp$apply.temp..1..median.>=0.5,"LC","NOR")
  temp1<-temp[temp$Predicted=="LC",]
  result[i,"RF_dis_specificity"]<-mean(temp1$Predicted==temp1$Actual)
  temp2<-length(which(temp1$Actual=="LC"))
  result[i,"RF_dis_sensitivity"]<-temp2/139
  
  temp<-prob_validation_RF_list[[i]]
  temp<-data.frame(apply(temp, 1, median))
  rownames(temp)<-rownames(validation_x)
  temp$Actual<-NA
  temp$Actual[grep("LC",rownames(temp))]<-"LC"
  temp$Actual[grep("NOR",rownames(temp))]<-"NOR"
  temp$Predicted<-ifelse(temp$apply.temp..1..median.>=0.5,"LC","NOR")
  temp1<-temp[temp$Predicted=="LC",]
  result[i,"RF_vali_specificity"]<-mean(temp1$Predicted==temp1$Actual)
  temp2<-length(which(temp1$Actual=="LC"))
  result[i,"RF_vali_sensitivity"]<-temp2/26
  
  temp<-prob_discovery_SVM_list[[i]]
  temp<-data.frame(apply(temp, 1, median))
  rownames(temp)<-rownames(discovery_x)
  temp$Actual<-NA
  temp$Actual[grep("LC",rownames(temp))]<-"LC"
  temp$Actual[grep("NOR",rownames(temp))]<-"NOR"
  temp$Predicted<-ifelse(temp$apply.temp..1..median.>=0.5,"LC","NOR")
  temp1<-temp[temp$Predicted=="LC",]
  result[i,"SVM_dis_specificity"]<-mean(temp1$Predicted==temp1$Actual)
  temp2<-length(which(temp1$Actual=="LC"))
  result[i,"SVM_dis_sensitivity"]<-temp2/139
  
  temp<-prob_validation_SVM_list[[i]]
  temp<-data.frame(apply(temp, 1, median))
  rownames(temp)<-rownames(validation_x)
  temp$Actual<-NA
  temp$Actual[grep("LC",rownames(temp))]<-"LC"
  temp$Actual[grep("NOR",rownames(temp))]<-"NOR"
  temp$Predicted<-ifelse(temp$apply.temp..1..median.>=0.5,"LC","NOR")
  temp1<-temp[temp$Predicted=="LC",]
  result[i,"SVM_vali_specificity"]<-mean(temp1$Predicted==temp1$Actual)
  temp2<-length(which(temp1$Actual=="LC"))
  result[i,"SVM_vali_sensitivity"]<-temp2/26
}
result<-data.frame(result)
result$LR_mean<-apply(result[,c(1:4)],1,mean)  
result$RF_mean<-apply(result[,c(5:8)],1,mean)  
result$SVM_mean<-apply(result[,c(9:12)],1,mean)  
result$mean<-apply(result[,c(1:12)],1,mean)  
write.csv(result,"Combination Sensitivity and specificity summary.csv")
