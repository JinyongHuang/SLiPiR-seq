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

####RNA name----
RNA_name <- list.files(path = "../../../Normalization_RPM_colsum/", pattern = "cpm.txt$", full.names = TRUE, recursive = TRUE) 
RNA_name <- RNA_name[grep("cfRNA",RNA_name, invert = TRUE)]
RNA_name <- lapply(RNA_name, read.table, header=TRUE, sep="\t",row.names=1)
RNA_name <- lapply(RNA_name, row.names)
names(RNA_name)<- c("lsRNA",  "miRNA",  "msRNA",  "piRNA",  "rsRNA",  "snoRNA", "snRNA",  "tsRNA",  "ysRNA")

####Sample ID----
meta<-read.csv("../../../Metadata/Meta_with_clinical_and_Summary_good_samples.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="NOR_SZBA" |  sample_type=="LC_SZBU" | sample_type=="NOR_SZDW" |
                 sample_type=="BRC_NXYK" | sample_type=="CRC_NXYK" | sample_type=="GC_NXYK" | sample_type=="HCC_NXYK")
meta$group<-ifelse(meta$sample_type=="NOR_SZBA" | meta$sample_type=="NOR_SZDW","Controls","Cases")
meta$group<-factor(meta$group)
table(meta$sample_type)
table(meta$group)

###Discovery cohort
dis_case<-subset(meta,sample_type=="LC_SZDE" | sample_type=="BRC_NXYK" | sample_type=="CRC_NXYK" | sample_type=="GC_NXYK" | sample_type=="HCC_NXYK")
dis_control<-subset(meta,sample_type=="NOR_SZBA")
dis_ID<-c(dis_case$ID,dis_control$ID)
dis_type<-c(dis_case$group,dis_control$group)
ndiscovery<-length(dis_ID)
###Validation cohort
vali_case<-subset(meta,sample_type=="LC_SZBU")
vali_control<-subset(meta,sample_type=="NOR_SZDW")
vali_ID<-c(vali_case$ID,vali_control$ID)
vali_type<-c(vali_case$group,vali_control$group)
nvalidation<-length(vali_ID)

####Data for machine learning----
AllRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[,meta$ID]

####Significant RNA
sig_cfRNA<-read.table("../DE_cancers_vs._cancerfree/DESeq2 candidate results-cfRNA-Cases vs. Controls.txt", header = T)
sig_list<-list()
Distribution<-double()
for (RNA in names(RNA_name)) {
  sig_list[[RNA]]<-rownames(sig_cfRNA)[(rownames(sig_cfRNA) %in% RNA_name[[RNA]])]
  Distribution[RNA]=table(rownames(sig_cfRNA) %in% RNA_name[[RNA]])["TRUE"]
}
sum(Distribution)
sig_list[["goiRNA"]]<-c(sig_list[["msRNA"]],sig_list[["miRNA"]],sig_list[["snRNA"]],sig_list[["snoRNA"]],sig_list[["tsRNA"]])

####Filter by LASSO logistic regression----
goi=sig_list[["goiRNA"]]
subset<-AllRNA[goi,]
subset<-t(subset)
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
seeds<-1:100 #repeat 100 times
registerDoParallel(cl<-makeCluster(10))
coef<-foreach(seed=seeds, .combine="cbind",.packages=c("tidyverse","caret","glmnet")) %dopar% {
  set.seed(seed)
  print(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  ##Fit
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 1, nfolds = 10)
  LASSO <- glmnet(training_x, training_y, family = "binomial", alpha = 1, lambda = cvfit$lambda.1se)
  temp <- data.frame(matrix(coef(LASSO)))
  temp
}
stopCluster(cl)
rownames(coef)<-c("Intercept",colnames(subset))
colnames(coef)<-1:100
coef_filter <-data.frame(apply(coef[-1,], 1, function(x) length(which(x!=0))))
coef_filter$name<-rownames(coef_filter)
coef_filter <-coef_filter[coef_filter$apply.coef..1.....1..function.x..length.which.x....0... > 0,]
goi<-rownames(coef_filter)

ggplot(coef_filter, aes(x=apply.coef..1.....1..function.x..length.which.x....0...)) + 
  geom_histogram(color="black", fill="gray") + labs(x="Frequency of non-zero coefficient",y="cfRNA of interest",title="Cancer vs. cancer-free")+
  theme_classic()+theme(axis.text = element_text(colour = "black",size = 10),plot.title = element_text(hjust=0.5,face = "bold",size = 14), 
                        axis.title = element_text(size = 10), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))+
  annotate(geom="text", x=50, y=50, label="207 goi cfRNA",color="black",size=4)
ggsave(paste0("Histogram for LASSO coef distribution.pdf"), width = 4, height = 6)
ggsave(paste0("Histogram for LASSO coef distribution.png"), width = 4, height = 6, bg="white")

####distribution----
goi_list<-list()
for (RNA in names(RNA_name)) {
  goi_list[[RNA]]<-goi[(goi %in% RNA_name[[RNA]])]
}
goi_list[["cfRNA"]]<-c(goi_list[["msRNA"]],goi_list[["miRNA"]],goi_list[["snRNA"]],goi_list[["snoRNA"]],goi_list[["tsRNA"]])
goi_list<-Filter(function(x) length(x) > 0, goi_list)
orderd_name<-c("cfRNA","msRNA","miRNA","snRNA","snoRNA","tsRNA")
goi_list<- goi_list[match(orderd_name, names(goi_list))]

####Machine learning----
ML<-list()
seeds<-1:100 #repeat 100 times
for (RNA in names(goi_list)) {
  print(paste("Start working on",RNA))
  goi=goi_list[[RNA]]
  if (length(goi)<=1) next 
  subset<-AllRNA[goi,]
  subset<-t(subset)
  discovery_x <- subset[dis_ID,]
  discovery_y <- dis_type
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
    c(ACC_training,ACC_testing,AUC_training,AUC_testing,SenAt1Spe_training,SenAt1Spe_testing)
  }
  stopCluster(cl)
  colnames(result_seed)= c('ACC_training','ACC_testing','AUC_training','AUC_testing','SenAt1Spe_training','SenAt1Spe_testing')
  LR_df<-data.frame(result_seed)
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
    c(ACC_training,ACC_testing,AUC_training,AUC_testing,SenAt1Spe_training,SenAt1Spe_testing)
  }
  stopCluster(cl)
  colnames(result_seed)= c('ACC_training','ACC_testing','AUC_training','AUC_testing','SenAt1Spe_training','SenAt1Spe_testing')
  RF_df<-data.frame(result_seed)
  
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
    c(ACC_training,ACC_testing,AUC_training,AUC_testing,SenAt1Spe_training,SenAt1Spe_testing)
  }
  stopCluster(cl)
  colnames(result_seed)= c('ACC_training','ACC_testing','AUC_training','AUC_testing','SenAt1Spe_training','SenAt1Spe_testing')
  SVM_df<-data.frame(result_seed)
  
  results<-rbind(LR_df,RF_df,SVM_df)
  results$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  ML[[RNA]]<-results
}

####plot----
saveRDS(ML,"ML_results_Cancers vs. Cancer-free.rds")
AUC_Testing<-lapply(ML, "[[", "AUC_testing")
AUC_Testing<-data.frame(AUC_Testing)
AUC_Testing$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
AUC_Testing_median<-aggregate(AUC_Testing[,1:length(ML)],by=list(AUC_Testing$method),median)
AUC_Testing_df<-melt(AUC_Testing)
AUC_Testing_df$variable<-factor(AUC_Testing_df$variable,levels = orderd_name)
nGOI<- as.integer(lapply(goi_list, length))
names(nGOI) <- orderd_name
nGOI<- nGOI[nGOI>=2]
x_axix_text<-paste0(names(nGOI),"\n",nGOI)

ggplot(AUC_Testing_df, aes(x=variable, y= value))+geom_boxplot(outlier.shape = NA)+  geom_jitter(aes(color=method),shape=16, position=position_jitter(0.2),alpha=0.8) + 
  labs(x="", y = "AUC in 100 iterations",title="Cancers vs. Cancer-free")+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+
  scale_x_discrete(labels = x_axix_text)+
  theme(axis.text.y = element_text(colour = "black",size = 12), axis.text.x = element_text(colour = "black",size = 12), 
        legend.text = element_text(size = 12), legend.title = element_blank(), legend.position = c(0.1,0.12), legend.background = element_blank(),
        axis.title = element_text(size = 12), panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        plot.title = element_text(hjust=0.5,face="bold",size = 12))
ggsave(paste0("AUC boxplot of different cfRNA-Lung cancer.pdf"), width = 5, height = 5)
ggsave(paste0("AUC boxplot of different cfRNA-Lung cancer.png"), width = 5, height = 5, bg="white")

##Risk score----
seeds<-1:100
dis_type<-factor(dis_type,levels=c("Controls","Cases"))
goi=goi_list[["cfRNA"]]
subset<-AllRNA[goi,]
subset<-t(subset)
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type####LR Risk score
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
  LR <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
  prob_discovery <- predict(LR, discovery_x, type = "response")
  prob_discovery
}
stopCluster(cl)
####boxplot discovery
df_discovery<-data.frame(id=rownames(discovery_x),prob_result,row.names = rownames(discovery_x))
df_discovery_median<-apply(df_discovery[,-1], 1, median)
df_discovery_median<-sort(df_discovery_median)
df_discovery_plot<-melt(df_discovery,"id")
df_discovery_plot$id<-factor(df_discovery_plot$id,levels = names(df_discovery_median))
df_discovery_plot$actual=substr(df_discovery_plot$id,1,3)
df_discovery_plot$actual=gsub("_","",df_discovery_plot$actual)
df_discovery_plot$actual=gsub("LC","Lung cancer (N=139)",df_discovery_plot$actual)
df_discovery_plot$actual=gsub("NOR","Cancer-free (N=106)",df_discovery_plot$actual)
df_discovery_plot$actual=gsub("BRC","Breast cancer (N=30)",df_discovery_plot$actual)
df_discovery_plot$actual=gsub("CRC","Colorectal cancer (N=37)",df_discovery_plot$actual)
df_discovery_plot$actual=gsub("GC","Gastric cancer (N=55)",df_discovery_plot$actual)
df_discovery_plot$actual=gsub("HCC","Liver cancer (N=15)",df_discovery_plot$actual)
df_discovery_plot$actual=factor(df_discovery_plot$actual,levels = c("Cancer-free (N=106)","Breast cancer (N=30)","Colorectal cancer (N=37)",
                                                                    "Gastric cancer (N=55)","Liver cancer (N=15)","Lung cancer (N=139)"))
discovery_plot<-ggplot(df_discovery_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.1)+ scale_colour_manual(values=c("#4E62AB","#469EB4","#87CFA4","#F0D43A","#F57547","#D6404E"))+
  labs(x=paste0("Study subjects (N=",length(discovery_y),")"), y="Cancer risk score")+ geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic()+theme(axis.text.x = element_blank(), axis.text.y = element_text(colour = "black"), text = element_text(size = 12), 
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=1), 
                        legend.position = "right", legend.title = element_blank(),legend.background = element_blank())
discovery_plot
ggsave(paste0("Cancers vs. Cancer-free-discovery-logistic regression.pdf"), width = 16, height = 3)
ggsave(paste0("Cancers vs. Cancer-free-discovery-logistic regression.png"), width = 16, height = 3, bg="white")

##median risk score boxplot
df_median<-as.data.frame(df_discovery_median)
df_median$actual=substr(rownames(df_median),1,3)
df_median$actual=gsub("_","",df_median$actual)
df_median$actual=factor(df_median$actual,levels = c("NOR","BRC","CRC","GC","HCC","LC"))

plot<-ggplot(df_median, aes(x=actual, y=df_discovery_median)) + geom_boxplot(outlier.shape=NA)+ geom_jitter(aes(colour = actual),shape=16, position=position_jitter(0.2),size=0.6)+
  scale_colour_manual(values=c("#4E62AB","#469EB4","#87CFA4","#F0D43A","#F57547","#D6404E"))+labs(x="", y="Cancer risk score")+
  theme_classic()+theme(text = element_text(size = 8), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black",angle = 30,hjust=1),
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.15), legend.position = "none",
                        axis.line = element_line(linewidth=0.15))
ggsave(paste0("Cancers median risk score boxplot by group.pdf"), width = 4, height = 4, units="cm")
ggsave(paste0("Cancers median risk score boxplot by group.png"), width = 4, height = 4, units="cm", bg="white")

####Different cancers goiRNA----
BRC<-readRDS("../ML_BRC_vs._others/ML_results_Breast cancer vs. Non-breast cancer.rds")
CRC<-readRDS("../ML_CRC_vs._others/ML_results_Colorectal cancer vs. Non-colorectal cancer.rds")
GC<-readRDS("../ML_GC_vs._others/ML_results_Gastric cancer vs. Non-gastric cancer.rds")
HCC<-readRDS("../ML_HCC_vs._others/ML_results_Liver cancer vs. Non-liver cancer.rds")
LC<-readRDS("../ML_LC_vs._others/ML_results_Lung cancer vs. Non-lung cancer.rds")
Different_cancers<-list(Cancers=ML,BRC=BRC,CRC=CRC,GC=GC,HCC=HCC,LC=LC)
Different_cancers<-lapply(Different_cancers, "[[", "cfRNA")
AUC_Testing<-lapply(Different_cancers, "[[", "AUC_testing")
AUC_Testing<-data.frame(AUC_Testing)
AUC_Testing$method<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
AUC_Testing_df<-melt(AUC_Testing)
names<-c("Cancers","BRC","CRC","GC","HCC","LC")
AUC_Testing_df$variable<-factor(AUC_Testing_df$variable,levels = names)
ggplot(AUC_Testing_df, aes(x=variable, y= value))+geom_boxplot(outlier.shape = NA)+  geom_jitter(aes(color=method),shape=16, position=position_jitter(0.2),alpha=0.8) + 
  labs(x="", y = "AUC in 100 iterations")+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+
  theme(axis.text.y = element_text(colour = "black",size = 12), axis.text.x = element_text(colour = "black",size = 12), 
        legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = c(0.12,0.12), legend.background = element_blank(),
        axis.title = element_text(size = 16), panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(paste0("AUC boxplot of different cancer types.pdf"), width = 4, height = 5)
ggsave(paste0("AUC boxplot of different cancer types.png"), width = 4, height = 5, bg="white")

##Risk score validation----
seeds<-1:100
dis_type<-factor(dis_type,levels=c("Controls","Cases"))
vali_type<-factor(vali_type,levels=c("Controls","Cases"))
goi=goi_list[["cfRNA"]]
subset<-AllRNA[goi,]
subset<-t(subset)
discovery_x <- subset[dis_ID,]
discovery_y <- dis_type
validation_x <- subset[vali_ID,]
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
  LR <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
  prob_validation <- predict(LR, validation_x, type = "response")
  prob_validation
}
stopCluster(cl)
####boxplot validation
df_validation<-data.frame(id=rownames(validation_x),prob_result,row.names = rownames(validation_x))
df_validation_median<-apply(df_validation[,-1], 1, median)
df_validation_median<-sort(df_validation_median)
saveRDS(df_validation_median,"Validation_Cancer_risk_score.rds")
df_validation_plot<-melt(df_validation,"id")
df_validation_plot$id<-factor(df_validation_plot$id,levels = names(df_validation_median))
df_validation_plot$actual=substr(df_validation_plot$id,1,3)
df_validation_plot$actual=gsub("_","",df_validation_plot$actual)
df_validation_plot$actual=gsub("LC","Lung cancer (N=26)",df_validation_plot$actual)
df_validation_plot$actual=gsub("NOR","Cancer-free (N=27)",df_validation_plot$actual)
df_validation_plot$actual=factor(df_validation_plot$actual,levels = c("Cancer-free (N=27)","Lung cancer (N=26)"))
validation_plot<-ggplot(df_validation_plot, aes(x=id, y=value,color=actual)) + geom_boxplot(outlier.size = 0.1)+ scale_colour_manual(values=c("#4E62AB","#D6404E"))+
  labs(x=paste0("Study subjects (N=",length(validation_y),")"), y="Cancer risk score")+ geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic()+theme(axis.text.x = element_blank(), axis.text.y = element_text(colour = "black"), text = element_text(size = 14),
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                        legend.position = c(0.18, 0.78), legend.title = element_blank(),legend.background = element_blank())
validation_plot