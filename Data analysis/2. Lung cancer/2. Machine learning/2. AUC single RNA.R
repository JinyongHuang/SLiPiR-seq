library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
#colorRampPalette(brewer.pal(9, "Blues"))(100)

TopN<-readRDS("../TopN/TopN_results.rds")
Boruta<-readRDS("../Boruta_10/Boruta_results.rds")
LASSO<-readRDS("../LASSO_10/LASSO_results.rds")

mRNA<-list(TopN[["mRNA"]],Boruta[["mRNA"]],LASSO[["mRNA"]])
lncRNA<-list(TopN[["lncRNA"]],Boruta[["lncRNA"]],LASSO[["lncRNA"]])
miRNA<-list(TopN[["miRNA"]],Boruta[["miRNA"]],LASSO[["miRNA"]])
piRNA<-list(TopN[["piRNA"]],Boruta[["piRNA"]],LASSO[["piRNA"]])
snRNA<-list(TopN[["snRNA"]],Boruta[["snRNA"]],LASSO[["snRNA"]])
snoRNA<-list(TopN[["snoRNA"]],Boruta[["snoRNA"]],LASSO[["snoRNA"]])
tsRNA<-list(TopN[["tsRNA"]],Boruta[["tsRNA"]],LASSO[["tsRNA"]])
rsRNA<-list(TopN[["rsRNA"]],Boruta[["rsRNA"]],LASSO[["rsRNA"]])
ysRNA<-list(TopN[["ysRNA"]],Boruta[["ysRNA"]],LASSO[["ysRNA"]])
RNA_list<-list(mRNA,lncRNA,miRNA,piRNA,snRNA,snoRNA,tsRNA,rsRNA,ysRNA)
names(RNA_list)<-c("mRNA","lncRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA","rsRNA","ysRNA")
RNA_list <- lapply(RNA_list, function(x) {names(x) <- c("TopN","Boruta","LASSO"); x })
RNA_result <- do.call(rbind, lapply(RNA_list, data.frame, stringsAsFactors = FALSE))

mRNA<-RNA_result[grep("mRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
lncRNA<-RNA_result[grep("lncRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
miRNA<-RNA_result[grep("miRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
piRNA<-RNA_result[grep("piRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
snRNA<-RNA_result[grep("snRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
snoRNA<-RNA_result[grep("snoRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
tsRNA<-RNA_result[grep("tsRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
rsRNA<-RNA_result[grep("rsRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
ysRNA<-RNA_result[grep("ysRNA",rownames(RNA_result)),grep("AUC",colnames(RNA_result))]
AUC_result<-list(mRNA,lncRNA,miRNA,piRNA,snRNA,snoRNA,tsRNA,rsRNA,ysRNA)
names(AUC_result)<-names(RNA_list)


####9 cfRNA good test and validation
mRNA_good<-mRNA[,grep("LASSO",colnames(mRNA))]
colnames(mRNA_good)<-c("mRNA_training","mRNA_test","mRNA_validation")
lncRNA_good<-lncRNA[,grep("LASSO",colnames(lncRNA))]
colnames(lncRNA_good)<-c("lncRNA_training","lncRNA_test","lncRNA_validation")
miRNA_good<-miRNA[,grep("LASSO",colnames(miRNA))]
colnames(miRNA_good)<-c("miRNA_training","miRNA_test","miRNA_validation")
piRNA_good<-piRNA[,grep("LASSO",colnames(piRNA))]
colnames(piRNA_good)<-c("piRNA_training","piRNA_test","piRNA_validation")
snRNA_good<-snRNA[,grep("LASSO",colnames(snRNA))]
colnames(snRNA_good)<-c("snRNA_training","snRNA_test","snRNA_validation")
snoRNA_good<-snoRNA[,grep("LASSO",colnames(snoRNA))]
colnames(snoRNA_good)<-c("snoRNA_training","snoRNA_test","snoRNA_validation")
tsRNA_good<-tsRNA[,grep("LASSO",colnames(tsRNA))]
colnames(tsRNA_good)<-c("tsRNA_training","tsRNA_test","tsRNA_validation")
rsRNA_good<-rsRNA[,grep("LASSO",colnames(rsRNA))]
colnames(rsRNA_good)<-c("rsRNA_training","rsRNA_test","rsRNA_validation")
ysRNA_good<-ysRNA[,grep("LASSO",colnames(ysRNA))]
colnames(ysRNA_good)<-c("ysRNA_training","ysRNA_test","ysRNA_validation")
allRNA<-cbind(mRNA_good,lncRNA_good,miRNA_good,piRNA_good,snRNA_good,snoRNA_good,tsRNA_good,rsRNA_good,ysRNA_good)
allRNA$ML<-c(rep("LR",100),rep("RF",100),rep("SVM",100))

df<-melt(allRNA,"ML")
df$variable<-gsub("_","\n",df$variable)
df$RNA<-c(rep("mRNA",900),rep("lncRNA",900),rep("miRNA",900),rep("piRNA",900),rep("snRNA",900),rep("snoRNA",900),rep("tsRNA",900),rep("rsRNA",900),rep("ysRNA",900))
df$RNA<-factor(df$RNA,levels = c("mRNA","lncRNA","miRNA","piRNA","snRNA","snoRNA","tsRNA","rsRNA","ysRNA"))

df$Group<-c(rep("Training",300),rep("Test",300),rep("Validation",300))
df$Group<-factor(df$Group,levels = c("Training","Test","Validation"))

ggplot(df, aes(x=Group, y= value))+geom_boxplot(aes(fill=Group),outlier.shape = NA, lwd=0.3,fatten=0.8)+  geom_jitter(aes(color=ML),shape=16, position=position_jitter(0.2),alpha=0.8,size=0.2) + facet_wrap(~RNA,nrow=2,ncol=5)+
  labs(x="", y = "AUC in 100 iterations")+ theme_classic()+scale_colour_manual(values=c("#D6404E","#87CFA4","#4E62AB"))+scale_fill_manual(values=c("#636363","#BDBDBD","#F7F7F7"))+
  theme(axis.text.y = element_text(colour = "black",size = 6), axis.text.x = element_blank(),  
        strip.text.x = element_text(colour = "black",size = 6), strip.background = element_blank(),axis.ticks.x = element_blank(),
        legend.text = element_text(size = 6), legend.title = element_blank(), legend.position = c(0.91,0.4), legend.background = element_blank(),
        axis.title = element_text(size = 6), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth=0.15))
ggsave(paste0("AUC boxplot-9 cfRNA good features.pdf"), width = 12, height = 6, units="cm")
ggsave(paste0("AUC boxplot-9 cfRNA good features.png"), width = 12, height = 6, units="cm", bg="white")

 
