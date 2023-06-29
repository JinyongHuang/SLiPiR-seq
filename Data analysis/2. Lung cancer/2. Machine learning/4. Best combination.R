library(tidyverse)
library(caret)#Sample partition, Support Vector Machine 
library(glmnet)#LASSO logistic regression
library(randomForest)#RF
library(MASS)#LDA
library(ROCR)#ROC
library(pROC)#ROC
library(reshape2)
library(ggplot2)
library(ggpubr)
library(Rtsne)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(RColorBrewer)
#colorRampPalette(brewer.pal(9, "Blues"))(100)

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
goi_snoRNA<-readLines("../TopN/goi/snoRNA Top N selected name.txt")
goi_tsRNA<-readLines("../LASSO/goi/20/tsRNA name.txt")

goi<-c(goi_msRNA,goi_miRNA,goi_snoRNA,goi_tsRNA)
selected="ms+mi+sno+ts"
####Data for machine learning----
AllRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[goi,meta$ID]

####best combination----
seeds<-1:100 #repeat 100 times
top_rpm <- t(AllRNA)
discovery_x <- top_rpm[dis_ID,]
discovery_y <- dis_type
validation_x <- top_rpm[vali_ID,]
validation_y <- vali_type

###Logistic regression
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
  c(ACC_training,ACC_testing,ACC_validation,
    AUC_training,AUC_testing,AUC_validation,
    SenAt1Spe_training,SenAt1Spe_testing,SenAt1Spe_validation)}
stopCluster(cl)
colnames(result_seed)= c('ACC_training','ACC_testing','ACC_validation',
                         'AUC_training','AUC_testing','AUC_validation',
                         'SenAt1Spe_training','SenAt1Spe_testing','SenAt1Spe_validation')
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
  c(ACC_training,ACC_testing,ACC_validation,
    AUC_training,AUC_testing,AUC_validation,
    SenAt1Spe_training,SenAt1Spe_testing,SenAt1Spe_validation)}
stopCluster(cl)
colnames(result_seed)= c('ACC_training','ACC_testing','ACC_validation',
                         'AUC_training','AUC_testing','AUC_validation',
                         'SenAt1Spe_training','SenAt1Spe_testing','SenAt1Spe_validation')
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
  c(ACC_training,ACC_testing,ACC_validation,
    AUC_training,AUC_testing,AUC_validation,
    SenAt1Spe_training,SenAt1Spe_testing,SenAt1Spe_validation)
}
stopCluster(cl)
colnames(result_seed)= c('ACC_training','ACC_testing','ACC_validation',
                         'AUC_training','AUC_testing','AUC_validation',
                         'SenAt1Spe_training','SenAt1Spe_testing','SenAt1Spe_validation')
SVM_df<-result_seed

Best<-data.frame(rbind(LR_df,RF_df,SVM_df))

summary(LR_df[,"AUC_testing"])
summary(LR_df[,"AUC_validation"])
# summary(RF_df[,"AUC_testing"])
# summary(RF_df[,"AUC_validation"])
# summary(SVM_df[,"AUC_testing"])
# summary(SVM_df[,"AUC_validation"])

df<-Best[,c(5,6,10)]
df<-melt(df,"method")
df$variable<-factor(c(rep("Testing",300),rep("Validation",300)),levels = c("Testing","Validation"))
ggplot(df, aes(x=variable, y= value,colour=variable))+geom_boxplot(outlier.shape = NA)+  geom_jitter(shape=16, position=position_jitter(0.2)) + facet_wrap(~method, scale="free")+
  labs(x="", y = "AUC in 100 iterations")+ theme_classic()+scale_colour_brewer(palette="Set1")+
  theme(axis.text.y = element_text(colour = "black",size = 14), axis.text.x = element_text(colour = "black",size = 14),
        strip.text.x = element_text(colour = "black",size = 14,face = "bold"), legend.position = "none", 
        axis.title = element_text(size = 16), plot.title = element_text(size = 16,hjust = 0.5,face = "bold"), panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave(paste0(selected," cfRNA panel - AUC.pdf"), width = 7, height = 7)
ggsave(paste0(selected," cfRNA panel - AUC.png"), width = 7, height = 7, bg="white")

#####cancer stage ----
dis_case_early<-subset(meta,sample_type=="LC_SZDE" & Stage=="Early")
dis_case_late<-subset(meta,sample_type=="LC_SZDE" & Stage=="Late")
dis_ID_early<-dis_case_early$ID
dis_ID_late<-dis_case_late$ID
discovery_x_early <- top_rpm[dis_ID_early,]
discovery_x_late <- top_rpm[dis_ID_late,]
discovery_y_early <- dis_case_early$group
discovery_y_late <- dis_case_late$group
####STAGE risk score
registerDoParallel(cl<-makeCluster(10))
prob_result<-foreach(seed=seeds, .combine="cbind", .packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  ##Fit
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
  Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)
  ##Discovery cohort
  prob_dis_early <- predict(Ridge, discovery_x_early, type = "response")
  prob_dis_late <- predict(Ridge, discovery_x_late, type = "response")
  c(prob_dis_early,prob_dis_late)
}
stopCluster(cl)
prob_discovery_early<-prob_result[1:length(discovery_y_early),]
prob_discovery_late<-prob_result[(length(discovery_y_early)+1):(length(discovery_y_early)+length(discovery_y_late)),]
rownames(prob_discovery_early)=rownames(discovery_x_early)
rownames(prob_discovery_late)=rownames(discovery_x_late)
####risk score boxplot 
df_early<-data.frame(id=rownames(discovery_x_early),prob_discovery_early,row.names = rownames(discovery_x_early))
df_late<-data.frame(id=rownames(discovery_x_late),prob_discovery_late,row.names = rownames(discovery_x_late))
df<-rbind(df_early,df_late)
df_median<-apply(df[,-1], 1, median)
df_median<-sort(df_median)
df_plot<-melt(df,"id")
df_plot$id<-factor(df_plot$id,levels = names(df_median))
df_plot$AJCCstage<-rep(c(dis_case_early$AJCC.Stage,dis_case_late$AJCC.Stage),100)
df_plot$AJCCstage=str_c("Stage ",df_plot$AJCCstage)
df_plot$AJCCstage=factor(df_plot$AJCCstage,levels = c("Stage I","Stage II","Stage III","Stage IV"))
ggplot(df_plot, aes(x=id, y=value,color=AJCCstage)) + geom_boxplot(outlier.size = 0.1)+ scale_color_manual(values = c("#FB7B5B","#F5553C","#E32F27","#C2161B"))+
  labs(x=paste0("Lung cancer patients (N=139)"), y="Lung cancer risk scores", title=paste0(""))+ geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic()+theme(axis.text.x = element_blank(),  axis.text.y = element_text(color = "black"), text = element_text(size = 8), 
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.line = element_line(linewidth=0.15),
                        legend.position = c(0.88, 0.28), legend.title = element_blank(),legend.background = element_blank())
ggsave(paste0("Risk score of different AJCC stage lung cancer.pdf"), width = 5, height = 2.5)
ggsave(paste0("Risk score of different AJCC stage lung cancer.png"), width = 5, height = 2.5, bg="white")

"#FFF5F0" "#FEE2D5" "#FCC3AB" "#FC9F81" "#FB7B5B" "#F5553C" "#E32F27" "#C2161B" "#9E0D14" "#67000D"
####STAGE AUC
seeds<-1:100
registerDoParallel(cl<-makeCluster(10))
AUC_result<-foreach(seed=seeds, .combine="rbind", .packages=c("tidyverse","caret","glmnet","ROCR","pROC")) %dopar% {
  set.seed(seed)
  splitSample <- createDataPartition(discovery_y, p = 0.8, list = FALSE)
  training_x <- discovery_x[splitSample,]
  training_y <- discovery_y[splitSample]
  testing_x <- discovery_x[-splitSample,]
  testing_y <- discovery_y[-splitSample]
  ##Fit
  cvfit <- cv.glmnet(training_x, training_y, family = "binomial", alpha = 0, nfolds = 10)
  Ridge <- glmnet(training_x, training_y, family = "binomial", alpha = 0, lambda = cvfit$lambda.min)

  testing<-data.frame(ID=rownames(testing_x))
  testing<-merge(testing,meta[,c("ID","group","Stage","AJCC.Stage")],by="ID")
  testing_normal<-testing[testing$group=="Healthy",]
  
  testing_early<-testing[testing$Stage=="Early",]
  testing_early<-na.omit(testing_early)
  testing_x_early<-testing_x[c(testing_early$ID,testing_normal$ID),]
  testing_y_early<-c(testing_early$group,testing_normal$group)
  prob_early <- predict(Ridge, testing_x_early, type = "response")
  pred_early <- prediction(prob_early,testing_y_early)
  AUC_early <- round(attr(performance(pred_early, "auc"), "y.values")[[1]],3)
  
  testing_late<-testing[testing$Stage=="Late",]
  testing_late<-na.omit(testing_late)
  testing_x_late<-testing_x[c(testing_late$ID,testing_normal$ID),]
  testing_y_late<-c(testing_late$group,testing_normal$group)
  prob_late <- predict(Ridge, testing_x_late, type = "response")
  pred_late <- prediction(prob_late,testing_y_late)
  AUC_late <- round(attr(performance(pred_late, "auc"), "y.values")[[1]],3)
  c(AUC_early,AUC_late)
}
stopCluster(cl)
colnames(AUC_result)<-c("Early","Late")
AUC_result<-data.frame(AUC_result)
AUC_result$seed<-1:100
df_plot<-melt(AUC_result,"seed")
df_plot$variable<-gsub("Early","Early stage",df_plot$variable)
df_plot$variable<-gsub("Late","Late stage",df_plot$variable)
ggplot(df_plot, aes(x=variable, y=value,color=variable)) + geom_boxplot(outlier.shape = NA)+ geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  scale_color_manual(values=c("#D6404E","#9E0142"))+labs(x="", y="AUC in 100 repeats", title=paste0(""))+
  theme_classic()+theme(text = element_text(size = 8), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), 
                        legend.position = "none", axis.text = element_text(color = "black"),axis.line = element_line(linewidth=0.15))
ggsave("AUC boxplot beween early and late stage.png",device = "png",bg="white",width = 4, height = 6,units = "cm")
ggsave("AUC boxplot beween early and late stage.pdf",device = "pdf",width = 4, height = 6,units = "cm")
