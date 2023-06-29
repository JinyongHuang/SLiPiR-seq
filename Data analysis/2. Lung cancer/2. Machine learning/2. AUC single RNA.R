library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)

TopN<-readRDS("../TopN/TopN_results.rds")
Boruta<-readRDS("../Boruta/Boruta_results.rds")
LASSO<-readRDS("../LASSO/LASSO_results.rds")

msRNA<-list(TopN[["msRNA"]],Boruta[["msRNA"]][["20"]],LASSO[["msRNA"]][["20"]])
lsRNA<-list(TopN[["lsRNA"]],Boruta[["lsRNA"]][["20"]],LASSO[["lsRNA"]][["20"]])
miRNA<-list(TopN[["miRNA"]],Boruta[["miRNA"]][["20"]],LASSO[["miRNA"]][["20"]])
piRNA<-list(TopN[["piRNA"]],Boruta[["piRNA"]][["20"]],LASSO[["piRNA"]][["20"]])
snRNA<-list(TopN[["snRNA"]],Boruta[["snRNA"]][["20"]],LASSO[["snRNA"]][["20"]])
snoRNA<-list(TopN[["snoRNA"]],Boruta[["snoRNA"]][["20"]],LASSO[["snoRNA"]][["20"]])
tsRNA<-list(TopN[["tsRNA"]],Boruta[["tsRNA"]][["20"]],LASSO[["tsRNA"]][["20"]])
rsRNA<-list(TopN[["rsRNA"]],Boruta[["rsRNA"]][["20"]],LASSO[["rsRNA"]][["20"]])
ysRNA<-list(TopN[["ysRNA"]],Boruta[["ysRNA"]][["20"]],LASSO[["ysRNA"]][["20"]])
RNA_list<-list(msRNA,lsRNA,miRNA,piRNA,snRNA,snoRNA,tsRNA,rsRNA,ysRNA)
names(RNA_list)<-c("msRNA","lsRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA","rsRNA","ysRNA")
RNA_list <- lapply(RNA_list, function(x) {names(x) <- c("TopN","Boruta","LASSO"); x })
RNA_result <- do.call(rbind, lapply(RNA_list, data.frame, stringsAsFactors = FALSE))

msRNA<-RNA_result[grep("msRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
lsRNA<-RNA_result[grep("lsRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
miRNA<-RNA_result[grep("miRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
piRNA<-RNA_result[grep("piRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
snRNA<-RNA_result[grep("snRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
snoRNA<-RNA_result[grep("snoRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
tsRNA<-RNA_result[grep("tsRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
rsRNA<-RNA_result[grep("rsRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
ysRNA<-RNA_result[grep("ysRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
AUC_result<-list(msRNA,lsRNA,miRNA,piRNA,snRNA,snoRNA,tsRNA,rsRNA,ysRNA)
names(AUC_result)<-names(RNA_list)
AUC_result_df<-lapply(AUC_result, melt)

plot_list<-list()
for (RNA in names(AUC_result)) {
  subset<-AUC_result_df[[RNA]]
  subset$ML<-c(rep("LR",100),rep("RF",100),rep("SVM",100))
  subset$method<-matrix(unlist(str_split(subset$variable,"\\.")),ncol=nrow(subset), nrow=2)[1,]
  subset$method<-factor(subset$method,levels =c("TopN","Boruta","LASSO"))
  subset$variable<-gsub("\\.","\n",subset$variable)
  subset$variable<-gsub("AUC_training","Training",subset$variable)
  subset$variable<-gsub("AUC_testing","Testing",subset$variable)
  subset$variable<-gsub("AUC_validation","Validation",subset$variable)
  subset$variable<-factor(subset$variable,levels =c("TopN\nTraining","TopN\nTesting","TopN\nValidation","Boruta\nTraining","Boruta\nTesting","Boruta\nValidation","LASSO\nTraining","LASSO\nTesting","LASSO\nValidation"))
  plot<-ggplot(subset, aes(x=variable, y= value))+geom_boxplot(outlier.shape = NA)+  geom_jitter(aes(color=ML),shape=16, position=position_jitter(0.2),alpha=0.8) + 
    labs(x="", y = "AUC in 100 iterations",title=paste0(RNA))+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+
    geom_vline(xintercept=3.5, linetype="dashed", color = "black") + geom_vline(xintercept=6.5, linetype="dashed", color = "black") +
    theme(axis.text.y = element_text(colour = "black",size = 12), axis.text.x = element_text(colour = "black",size = 12), 
          legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = c(0.08,0.10), legend.background = element_blank(),
          axis.title = element_text(size = 16), panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          plot.title = element_text(hjust=0.5,face="bold"))
  ggsave(paste0("AUC boxplot-", RNA,".pdf"), width = 7, height = 7)
  ggsave(paste0("AUC boxplot-", RNA,".png"), width = 7, height = 7, bg="white")
  plot_list[[RNA]]<-plot
}

####9 cfRNA good test and validation
msRNA_good<-msRNA[,grep("LASSO",colnames(msRNA))]
colnames(msRNA_good)<-c("msRNA_training","msRNA_test","msRNA_validation")
lsRNA_good<-lsRNA[,grep("TopN",colnames(lsRNA))]
colnames(lsRNA_good)<-c("lsRNA_training","lsRNA_test","lsRNA_validation")
miRNA_good<-miRNA[,grep("Boruta",colnames(miRNA))]
colnames(miRNA_good)<-c("miRNA_training","miRNA_test","miRNA_validation")
piRNA_good<-piRNA[,grep("Boruta",colnames(piRNA))]
colnames(piRNA_good)<-c("piRNA_training","piRNA_test","piRNA_validation")
snRNA_good<-snRNA[,grep("LASSO",colnames(snRNA))]
colnames(snRNA_good)<-c("snRNA_training","snRNA_test","snRNA_validation")
snoRNA_good<-snoRNA[,grep("TopN",colnames(snoRNA))]
colnames(snoRNA_good)<-c("snoRNA_training","snoRNA_test","snoRNA_validation")
tsRNA_good<-tsRNA[,grep("LASSO",colnames(tsRNA))]
colnames(tsRNA_good)<-c("tsRNA_training","tsRNA_test","tsRNA_validation")
rsRNA_good<-rsRNA[,grep("LASSO",colnames(rsRNA))]
colnames(rsRNA_good)<-c("rsRNA_training","rsRNA_test","rsRNA_validation")
ysRNA_good<-ysRNA[,grep("TopN",colnames(ysRNA))]
colnames(ysRNA_good)<-c("ysRNA_training","ysRNA_test","ysRNA_validation")
allRNA<-cbind(msRNA_good,lsRNA_good,miRNA_good,piRNA_good,snRNA_good,snoRNA_good,tsRNA_good,rsRNA_good,ysRNA_good)
allRNA$ML<-c(rep("LR",100),rep("RF",100),rep("SVM",100))

df<-melt(allRNA,"ML")
df$variable<-gsub("_","\n",df$variable)
df$RNA<-c(rep("msRNA",900),rep("lsRNA",900),rep("miRNA",900),rep("piRNA",900),rep("snRNA",900),rep("snoRNA",900),rep("tsRNA",900),rep("rsRNA",900),rep("ysRNA",900))
#df$RNA<-factor(df$RNA,levels = c("msRNA","lsRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA","rsRNA","ysRNA"))
df$RNA<-gsub("msRNA","msRNA (LASSO)",df$RNA)
df$RNA<-gsub("lsRNA","lsRNA (TopN)",df$RNA)
df$RNA<-gsub("miRNA","miRNA (Boruta)",df$RNA)
df$RNA<-gsub("piRNA","piRNA (Boruta)",df$RNA)
df$RNA<-gsub("snRNA","snRNA (LASSO)",df$RNA)
df$RNA<-gsub("snoRNA","snoRNA (TopN)",df$RNA)
df$RNA<-gsub("tsRNA","tsRNA (LASSO)",df$RNA)
df$RNA<-gsub("rsRNA","rsRNA (LASSO)",df$RNA)
df$RNA<-gsub("ysRNA","ysRNA (TopN)",df$RNA)
df$RNA<-factor(df$RNA,levels = c("msRNA (LASSO)","lsRNA (TopN)","miRNA (Boruta)",
                                 "piRNA (Boruta)","snRNA (LASSO)","snoRNA (TopN)",
                                 "tsRNA (LASSO)","rsRNA (LASSO)","ysRNA (TopN)"))
df$Group<-c(rep("Training",300),rep("Test",300),rep("Validation",300))
df$Group<-factor(df$Group,levels = c("Training","Test","Validation"))

ggplot(df, aes(x=Group, y= value))+geom_boxplot(aes(fill=Group),outlier.shape = NA)+  geom_jitter(aes(color=ML),shape=16, position=position_jitter(0.2),alpha=0.8,size=0.2) + facet_wrap(~RNA,nrow=2,ncol=5)+
  labs(x="", y = "AUC in 100 iterations")+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+scale_fill_manual(values=c("#636363","#BDBDBD","#F7F7F7"))+
  theme(axis.text.y = element_text(colour = "black",size = 6), axis.text.x = element_blank(),  
        strip.text.x = element_text(colour = "black",size = 6), strip.background = element_blank(),
        legend.text = element_text(size = 6), legend.title = element_blank(), legend.position = c(0.91,0.25), legend.background = element_blank(),
        axis.title = element_text(size = 6), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth=0.15))
ggsave(paste0("AUC boxplot-9 cfRNA good features.pdf"), width = 12, height = 5, units="cm")
ggsave(paste0("AUC boxplot-9 cfRNA good features.png"), width = 12, height = 5, units="cm", bg="white")

 