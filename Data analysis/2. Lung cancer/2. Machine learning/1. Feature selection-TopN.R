library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(MASS)#LDA
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
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

####Significant RNA
RNA_name <- list.files(path = "../../../Normalization_RPM_colsum/", pattern = "cpm.txt$", full.names = TRUE, recursive = TRUE) 
RNA_name <- RNA_name[grep("cfRNA",RNA_name, invert = TRUE)]
RNA_name <- lapply(RNA_name, read.table, header=TRUE, sep="\t",row.names=1)
RNA_name <- lapply(RNA_name, row.names)
names(RNA_name)<- c("lncRNA",  "miRNA",  "mRNA",  "piRNA",  "rsRNA",  "snoRNA", "snRNA",  "tsRNA",  "ysRNA")
sig_cfRNA<-read.table("../../DE_LC_SZDE_vs_NOR_SZBA/DESeq2 candidate results-cfRNA-Lung cancer vs. Healthy.txt", header = T)
sig_cfRNA<-sig_cfRNA[sig_cfRNA$log2FoldChange>=0,]
sig_list<-list()
Distribution<-double()
for (RNA in names(RNA_name)) {
  sig_list[[RNA]]<-rownames(sig_cfRNA)[(rownames(sig_cfRNA) %in% RNA_name[[RNA]])]
  Distribution[RNA]=table(rownames(sig_cfRNA) %in% RNA_name[[RNA]])["TRUE"]
}
sum(Distribution)

####Data for machine learning----
cfRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
mRNA<-cfRNA[sig_list[["mRNA"]],meta$ID]
lncRNA<-cfRNA[sig_list[["lncRNA"]],meta$ID]
miRNA<-cfRNA[sig_list[["miRNA"]],meta$ID]
piRNA<-cfRNA[sig_list[["piRNA"]],meta$ID]
tsRNA<-cfRNA[sig_list[["tsRNA"]],meta$ID]
rsRNA<-cfRNA[sig_list[["rsRNA"]],meta$ID]
ysRNA<-cfRNA[sig_list[["ysRNA"]],meta$ID]
snRNA<-cfRNA[sig_list[["snRNA"]],meta$ID]
snoRNA<-cfRNA[sig_list[["snoRNA"]],meta$ID]
RNA_list<-list(mRNA=mRNA,lncRNA=lncRNA,miRNA=miRNA,piRNA=piRNA,tsRNA=tsRNA,snRNA=snRNA,snoRNA=snoRNA,rsRNA=rsRNA,ysRNA=ysRNA)

####Machine learning general
seeds<-1:100 #repeat 100 times
####Ridge logistic regression----
registerDoParallel(cl<-makeCluster(10))
LR_result<-list()
for (RNA in names(RNA_list)) {
  subset<-RNA_list[[RNA]]
  goi<-c(2:nrow(subset))
  if (nrow(subset)>1000) {goi<-seq(2,nrow(subset),10)}#rsRNA
  result_goi<-data.frame(AUC_training_mean=rep(0,length(goi)),AUC_testing_mean=rep(0,length(goi)),
                         AUC_training_sd=rep(0,length(goi)),AUC_testing_sd=rep(0,length(goi)),
                         row.names = goi)
  for (N in goi) {
    print(paste("Start working on top",N,RNA))
    top_N<-sig_list[[RNA]][1:N] #ranked by pvalue
    top_rpm<-subset[top_N,]
    top_rpm <- t(top_rpm)
    discovery_x <- top_rpm[dis_ID,]
    discovery_y <- dis_type
    result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","glmnet","ROCR")) %dopar% {
      set.seed(seed)
      splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
      training_x <- discovery_x[splitSample,]
      training_y <- discovery_y[splitSample]
      testing_x <- discovery_x[-splitSample,]
      testing_y <- discovery_y[-splitSample]
      ##Fit
      cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
      LR <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
      ##Training set
      prob <- predict(LR, training_x, type = "response")
      pred <- prediction(prob, training_y)
      AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
      ##Held-out testing set
      prob <- predict(LR, testing_x, type = "response")
      pred <- prediction(prob, testing_y)
      AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
      c(AUC_training,AUC_testing)
    }
    AUC_training_mean<-round(mean(result_seed[,1]),2)
    AUC_testing_mean<-round(mean(result_seed[,2]),2)
    AUC_training_sd<-round(sd(result_seed[,1]),2)
    AUC_testing_sd<-round(sd(result_seed[,2]),2)
    result_goi[as.character(N),]<-c(AUC_training_mean,AUC_testing_mean,AUC_training_sd,AUC_testing_sd)
  }
  LR_result[[RNA]]=result_goi
}
stopCluster(cl)

####Random forest----
registerDoParallel(cl<-makeCluster(10))
RF_result<-list()
for (RNA in names(RNA_list)) {
  subset<-RNA_list[[RNA]]
  goi<-c(2:nrow(subset))
  if (nrow(subset)>1000) {goi<-seq(2,nrow(subset),10)}#rsRNA
  result_goi<-data.frame(AUC_training_mean=rep(0,length(goi)),AUC_testing_mean=rep(0,length(goi)),
                         AUC_training_sd=rep(0,length(goi)),AUC_testing_sd=rep(0,length(goi)),
                         row.names = goi)
  for (N in goi) {
    print(paste("Start working on top",N,RNA))
    top_N<-sig_list[[RNA]][1:N] #ranked by pvalue
    top_rpm<-subset[top_N,]
    top_rpm <- t(top_rpm)
    discovery_x <- top_rpm[dis_ID,]
    discovery_y <- dis_type
    result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","randomForest","ROCR")) %dopar% {
      set.seed(seed)
      splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
      training_x <- discovery_x[splitSample,]
      training_y <- discovery_y[splitSample]
      testing_x <- discovery_x[-splitSample,]
      testing_y <- discovery_y[-splitSample]
      ##Fit
      RF = randomForest(x = training_x, y = training_y)
      ##Training set
      prob <- predict(RF, training_x, type = "prob")
      pred <- prediction(prob[,2], training_y)
      AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
      ##Held-out testing set
      prob <- predict(RF, testing_x, type = "prob")
      pred <- prediction(prob[,2], testing_y)
      AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
      c(AUC_training,AUC_testing)
    }
    AUC_training_mean<-round(mean(result_seed[,1]),2)
    AUC_testing_mean<-round(mean(result_seed[,2]),2)
    AUC_training_sd<-round(sd(result_seed[,1]),2)
    AUC_testing_sd<-round(sd(result_seed[,2]),2)
    result_goi[as.character(N),]<-c(AUC_training_mean,AUC_testing_mean,AUC_training_sd,AUC_testing_sd)
  }
  RF_result[[RNA]]=result_goi
}
stopCluster(cl)

####Support vector machine----
registerDoParallel(cl<-makeCluster(10))
SVM_result<-list()
for (RNA in names(RNA_list)) {
  subset<-RNA_list[[RNA]]
  goi<-c(2:nrow(subset))
  if (nrow(subset)>1000) {goi<-seq(2,nrow(subset),10)}#rsRNA
  result_goi<-data.frame(AUC_training_mean=rep(0,length(goi)),AUC_testing_mean=rep(0,length(goi)),
                         AUC_training_sd=rep(0,length(goi)),AUC_testing_sd=rep(0,length(goi)),
                         row.names = goi)
  for (N in goi) {
    print(paste("Start working on top",N,RNA))
    top_N<-sig_list[[RNA]][1:N] #ranked by pvalue
    top_rpm<-subset[top_N,]
    top_rpm <- t(top_rpm)
    discovery_x <- top_rpm[dis_ID,]
    discovery_y <- dis_type
    result_seed<-foreach(seed=seeds, .combine="rbind",.packages=c("tidyverse","caret","ROCR")) %dopar% {
      set.seed(seed)
      splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
      training_x <- discovery_x[splitSample,]
      training_y <- discovery_y[splitSample]
      testing_x <- discovery_x[-splitSample,]
      testing_y <- discovery_y[-splitSample]
      ##Fit
      trControl <- trainControl(method="cv", number=10, repeats=NA, p=0.8, classProbs = TRUE)# Set up Repeated k-fold Cross Validation
      SVM <- train(training_x, training_y, method="svmLinear", trControl=trControl, preProcess=c("center","scale"))
      ##Training set
      prob <- predict(SVM, training_x, type = "prob")
      pred <- prediction(prob[,2], training_y)
      AUC_training <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
      ##Held-out testing set
      prob <- predict(SVM, testing_x, type = "prob")
      pred <- prediction(prob[,2], testing_y)
      AUC_testing <- round(attr(performance(pred, "auc"), "y.values")[[1]],3)
      c(AUC_training,AUC_testing)
    }
    AUC_training_mean<-round(mean(result_seed[,1]),2)
    AUC_testing_mean<-round(mean(result_seed[,2]),2)
    AUC_training_sd<-round(sd(result_seed[,1]),2)
    AUC_testing_sd<-round(sd(result_seed[,2]),2)
    result_goi[as.character(N),]<-c(AUC_training_mean,AUC_testing_mean,AUC_training_sd,AUC_testing_sd)
  }
  SVM_result[[RNA]]=result_goi
}
stopCluster(cl)

####Identify the best top n----
dir.create(paste0(getwd(),"/Best"), showWarnings = FALSE)
AUC_list<-list()
topn<-double()
best_name<-list()
for (RNA in names(RNA_list)) {
  goi<-length(sig_list[[RNA]])
  topn_df<-data.frame(LR_training=LR_result[[RNA]]$AUC_training_mean,LR_testing=LR_result[[RNA]]$AUC_testing_mean,
                      RF_training=RF_result[[RNA]]$AUC_training_mean,RF_testing=RF_result[[RNA]]$AUC_testing_mean,
                      SVM_training=SVM_result[[RNA]]$AUC_training_mean,SVM_testing=SVM_result[[RNA]]$AUC_testing_mean,
                      topN=NA)
  if (RNA!="rsRNA") {topn_df$topN<-2:goi}
  if (RNA=="rsRNA") {topn_df$topN<-seq(2,goi,10)}#rsRNA
  
  topn_df$mean<-apply(topn_df[,-7],1,mean)
  topn_df$Training_mean<-apply(topn_df[,c(1,3,5)],1,mean)  
  topn_df$Testing_mean<-apply(topn_df[,c(2,4,6)],1,mean)  
  topn_df$LR_mean<-apply(topn_df[,1:2],1,mean)
  topn_df$RF_mean<-apply(topn_df[,3:4],1,mean)
  topn_df$SVM_mean<-apply(topn_df[,5:6],1,mean)
  
  topn_testing_df<-topn_df[order(-topn_df$Testing_mean),]
  print(paste0("Best mean testing AUC of ",RNA," is ",topn_testing_df$Testing_mean[1]," N=",topn_testing_df$topN[1]))
  write.csv(topn_testing_df,paste0("Best/",RNA," testing AUC.csv"))
  AUC_list[[RNA]]=topn_testing_df
  topn[[RNA]]=as.numeric(rownames(topn_testing_df)[1])
  Best_RNA<-sig_list[[RNA]][1:as.numeric(rownames(topn_testing_df)[1])]
  best_name[[RNA]]<-Best_RNA
  write_lines(Best_RNA,paste0("Best/Best_",RNA,"_name.txt"))
}

####Filter out from top N testing set AUC----
Filtered<-list()
for (RNA in names(RNA_list)) {
  temp<-AUC_list[[RNA]]
  temp<-temp[order(as.numeric(rownames(temp))),]
  temp<-temp[1:(topn[RNA]-1),]
  temp$effect<-"NA"
  for (i in 2:nrow(temp)){
    effect<-temp$Testing_mean[i]-temp$Testing_mean[i-1]
    if (effect<0) {temp$effect[i]="bad"}
    if (effect>0) {temp$effect[i]="good"}
  }
  bad=which(temp$effect=="bad")+1
  sig=sig_list[[RNA]][1:(topn[RNA])]
  if (length(bad)==0) {good=sig}
  if (length(bad)!=0) {good=sig[-bad]}
  Filtered[[RNA]]<-good
}
nFiltered=sapply(Filtered, function(x) length(x))

####All in one----
Best<-list()
seeds<-1:100 #repeat 100 times
for (RNA in names(RNA_list)[1:9]) {
  print(paste("Start working on",RNA))
  top_rpm<-RNA_list[[RNA]][Filtered[[RNA]],]
  top_rpm <- t(top_rpm)
  discovery_x <- top_rpm[dis_ID,]
  discovery_y <- dis_type
  validation_x <- top_rpm[vali_ID,]
  validation_y <- vali_type
  ###Ridge logistic regression
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
  
  Best[[RNA]]<-rbind(LR_df,RF_df,SVM_df)
}  
saveRDS(Best,"TopN_results.rds")

###plot
dir.create(paste0(getwd(),"/Filterout negative"), showWarnings = FALSE)
for (RNA in names(RNA_list)[1:9]) {
  n=length(Filtered[[RNA]])
  df<-data.frame(Best[[RNA]][,5:6])
  df$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  df<-melt(df)
  df$variable<-factor(c(rep("Testing",300),rep("Validation",300)),levels = c("Training","Testing","Validation"))
  ggplot(df, aes(x=variable, y= value,colour=variable))+geom_boxplot(outlier.shape = NA)+  geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(~method, scale="free")+
    scale_colour_brewer(palette="Set1")+labs(title=paste0(n, " ",RNA," selected by Top N") ,x="", y = "AUC in 100 iterations")+ theme_classic()+
    theme(axis.text.y = element_text(colour = "black",size = 14), axis.text.x = element_text(colour = "black",size = 14),
          strip.text.x = element_text(colour = "black",size = 14,face = "bold"), legend.position = "none", 
          axis.title = element_text(size = 16), plot.title = element_text(size = 16,hjust = 0.5,face = "bold"), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  ggsave(paste0("Filterout negative/Machine learning - Filtered ", n, " ",RNA," - AUC.pdf"), width = 7, height = 7)
  ggsave(paste0("Filterout negative/Machine learning - Filtered ", n, " ",RNA," - AUC.png"), width = 7, height = 7, bg="white")
}

###summary
dir.create(paste0(getwd(),"/goi"), showWarnings = FALSE)
for (RNA in names(RNA_list)) {
  goi=Filtered[[RNA]]
  writeLines(goi,paste0("goi/",RNA," Top N selected name.txt"))
  
  df<-data.frame(Best[[RNA]][,4:6])
  df$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  summary<-aggregate(df[,1:3], by=list(method=df$method), FUN = median)
  summary$RNA=RNA
  if (RNA==names(RNA_list)[1]){Summary=summary}
  if (RNA!=names(RNA_list)[1]){Summary=rbind(Summary,summary)}
}

publication<-matrix(ncol = 9, nrow = 10)
colnames(publication)<-names(RNA_list)
rownames(publication)<-c("LR_training","RF_training","SVM_training",
                         "LR_testing","RF_testing","SVM_testing",
                         "LR_validation","RF_validation","SVM_validation",
                         "nFeatures")
publication[1,]<-Summary$AUC_training[seq(1,nrow(Summary),3)]
publication[2,]<-Summary$AUC_training[seq(2,nrow(Summary),3)]
publication[3,]<-Summary$AUC_training[seq(3,nrow(Summary),3)]
publication[4,]<-Summary$AUC_testing[seq(1,nrow(Summary),3)]
publication[5,]<-Summary$AUC_testing[seq(2,nrow(Summary),3)]
publication[6,]<-Summary$AUC_testing[seq(3,nrow(Summary),3)]
publication[7,]<-Summary$AUC_validation[seq(1,nrow(Summary),3)]
publication[8,]<-Summary$AUC_validation[seq(2,nrow(Summary),3)]
publication[9,]<-Summary$AUC_validation[seq(3,nrow(Summary),3)]
publication[10,]<-nFiltered
publication<-publication[,c("mRNA",  "lncRNA",  "miRNA",  "piRNA",  
                            "snRNA",  "snoRNA", "tsRNA",  "rsRNA",  "ysRNA")]

write.csv(publication,"Top N selected AUC summary.csv")
