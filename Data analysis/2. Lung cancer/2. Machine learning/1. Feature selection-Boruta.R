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
names(RNA_name)<- c("lsRNA",  "miRNA",  "msRNA",  "piRNA",  "rsRNA",  "snoRNA", "snRNA",  "tsRNA",  "ysRNA")
sig_cfRNA<-read.table("../../DE_vsHealthy/DESeq2 significant results/LC_SE vs. NOR_SB/DESeq2 significant results-cfRNA-Lung cancer vs. Healthy.txt", header = T)
sig_cfRNA<-sig_cfRNA[sig_cfRNA$padj < 0.01 & sig_cfRNA$baseMean>=10 & sig_cfRNA$log2FoldChange>=1,]
sig_list<-list()
Distribution<-double()
for (RNA in names(RNA_name)) {
  sig_list[[RNA]]<-rownames(sig_cfRNA)[(rownames(sig_cfRNA) %in% RNA_name[[RNA]])]
  Distribution[RNA]=table(rownames(sig_cfRNA) %in% RNA_name[[RNA]])["TRUE"]
}
sum(Distribution)

####Data for machine learning----
cfRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
msRNA<-cfRNA[sig_list[["msRNA"]],meta$ID]
lsRNA<-cfRNA[sig_list[["lsRNA"]],meta$ID]
miRNA<-cfRNA[sig_list[["miRNA"]],meta$ID]
piRNA<-cfRNA[sig_list[["piRNA"]],meta$ID]
tsRNA<-cfRNA[sig_list[["tsRNA"]],meta$ID]
rsRNA<-cfRNA[sig_list[["rsRNA"]],meta$ID]
ysRNA<-cfRNA[sig_list[["ysRNA"]],meta$ID]
snRNA<-cfRNA[sig_list[["snRNA"]],meta$ID]
snoRNA<-cfRNA[sig_list[["snoRNA"]],meta$ID]
RNA_list<-list(msRNA=msRNA,lsRNA=lsRNA,miRNA=miRNA,piRNA=piRNA,tsRNA=tsRNA,snRNA=snRNA,snoRNA=snoRNA,rsRNA=rsRNA,ysRNA=ysRNA)

####Filter by Random forest Boruta----
goi_list<-list()
goi_N_list<-list()
condition<-seq(0,100,10)
for (RNA in names(RNA_list)) {
  print(paste("Start working on",RNA))
  subset<-t(RNA_list[[RNA]])
  discovery_x <- subset[dis_ID,]
  discovery_y <- dis_type
  seeds<-1:100 #repeat 100 times
  registerDoParallel(cl<-makeCluster(2))
  result_seed<-foreach(seed=seeds, .combine="c",.packages=c("tidyverse","caret","randomForest","Boruta")) %dopar% {
    set.seed(seed)
    print(seed)
    splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
    training_x <- discovery_x[splitSample,]
    training_y <- discovery_y[splitSample]
    testing_x <- discovery_x[-splitSample,]
    testing_y <- discovery_y[-splitSample]
    ##Fit
    boruta <- Boruta(training_x, training_y, doTrace = 2, maxRuns = 500)
    #plot(boruta, las = 2, cex.axis = 0.7)
    boruta <- TentativeRoughFix(boruta)
    goi<-names(boruta$finalDecision[boruta$finalDecision=="Confirmed"])
    goi
  }
  stopCluster(cl)
  for (N in condition) {
    top<-unique(result_seed)
    top<-top[table(result_seed) >= N]
    goi_N_list[[as.character(N)]]=top
    goi_N_list<-lapply(goi_N_list, function(x) gsub("\\.","-",x))
    goi_N_list<-lapply(goi_N_list, function(x) gsub("5-8S","5.8S",x))#5.8S rRNA
  }
  goi_list[[RNA]]<-goi_N_list
}
goi_list_summary<-data.frame(condition=seq(0,100,10))
goi_list_summary[names(RNA_list)]<-0
for (RNA in names(RNA_list)) {
  goi_list_summary[RNA]<-unlist(lapply(goi_list[[RNA]], function(x) length(x)))
}
dir.create(paste0(getwd(),"/goi"), showWarnings = FALSE)
for (N in condition) {
  dir.create(paste0(getwd(),"/goi/",N), showWarnings = FALSE)
  for (RNA in names(goi_list)) {
    goi=goi_list[[RNA]][[as.character(N)]]
    writeLines(goi,paste0("goi/",N,"/",RNA," name.txt"))
  }
}  
####All in one----
Best<-list()
condition<-20
seeds<-1:100 #repeat 100 times
for (RNA in names(goi_list)) {
  Best_N<-list()
    print(paste("Start working on",RNA))
    goi=goi_list[[RNA]][[as.character(condition)]]
    if (length(goi)<=1) break 
    subset<-RNA_list[[RNA]][goi,]
    subset<-t(subset)
    discovery_x <- subset[dis_ID,]
    discovery_y <- dis_type
    validation_x <- subset[vali_ID,]
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
    
    Best_N[[as.character(condition)]]<-rbind(LR_df,RF_df,SVM_df)
  
  Best[[RNA]]<-Best_N
  rm(Best_N)
}
saveRDS(Best,"Boruta_results.rds")

####summary----
con="20"
for (RNA in names(goi_list)) {
  df<-as.data.frame(Best[[RNA]])
  df$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  summary<-aggregate(df[,4:6], by=list(method=df$method), FUN = median)
  summary$RNA=RNA
  summary$features<-length(goi_list[[RNA]][[con]])
  if (RNA==names(goi_list)[1]){final=summary}
  if (RNA!=names(goi_list)[1]){final=rbind(final,summary)}
}

###publication
publication<-matrix(ncol = 9, nrow = 10)
colnames(publication)<-names(RNA_list)
rownames(publication)<-c("LR_training","RF_training","SVM_training",
                         "LR_testing","RF_testing","SVM_testing",
                         "LR_validation","RF_validation","SVM_validation",
                         "nFeatures")
publication[1,]<-final$X20.AUC_training[seq(1,nrow(final),3)]
publication[2,]<-final$X20.AUC_training[seq(2,nrow(final),3)]
publication[3,]<-final$X20.AUC_training[seq(3,nrow(final),3)]
publication[4,]<-final$X20.AUC_testing[seq(1,nrow(final),3)]
publication[5,]<-final$X20.AUC_testing[seq(2,nrow(final),3)]
publication[6,]<-final$X20.AUC_testing[seq(3,nrow(final),3)]
publication[7,]<-final$X20.AUC_validation[seq(1,nrow(final),3)]
publication[8,]<-final$X20.AUC_validation[seq(2,nrow(final),3)]
publication[9,]<-final$X20.AUC_validation[seq(3,nrow(final),3)]
nFeature=double()
rownames(goi_list_summary)<-goi_list_summary$condition
for (i in names(goi_list)) {
  nFeature[i]=goi_list_summary[con,i]
}
publication[10,]<-nFeature
publication<-publication[,c("msRNA",  "lsRNA",  "miRNA",  "piRNA",  
                            "snRNA",  "snoRNA", "tsRNA",  "rsRNA",  "ysRNA")]
write.csv(publication,"Boruta selected AUC summary.csv")

####plot
dir.create(paste0(getwd(),"/All in one"), showWarnings = FALSE)
for (RNA in names(RNA_list)) {
  n=length(goi_list[[RNA]][[con]])
  df<-data.frame(Best[[RNA]][[con]][,5:6])
  df$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  df<-melt(df)
  df$variable<-factor(c(rep("Testing",300),rep("Validation",300)),levels = c("Testing","Validation"))
  ggplot(df, aes(x=variable, y= value,colour=variable))+geom_boxplot(outlier.shape = NA)+  geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(~method, scale="free")+
    scale_colour_brewer(palette="Set1")+labs(title=paste0(n, " ",RNA," selected by Boruta") ,x="", y = "AUC in 100 iterations")+ theme_classic()+
    theme(axis.text.y = element_text(colour = "black",size = 14), axis.text.x = element_text(colour = "black",size = 14),
          strip.text.x = element_text(colour = "black",size = 14,face = "bold"), legend.position = "none", 
          axis.title = element_text(size = 16), plot.title = element_text(size = 16,hjust = 0.5,face = "bold"), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  ggsave(paste0("All in one/Machine learning - Boruta selected ", n, " ",RNA," - AUC.pdf"), width = 7, height = 7)
  ggsave(paste0("All in one/Machine learning - Boruta selected ", n, " ",RNA," - AUC.png"), width = 7, height = 7, bg="white")
}


####Median and IQR----
summary(Best[["miRNA"]][["20"]][1:100,"AUC_validation"])#Boruta selected miRNA in LR
summary(Best[["miRNA"]][["20"]][201:300,"AUC_validation"])#Boruta selected miRNA in SVM

