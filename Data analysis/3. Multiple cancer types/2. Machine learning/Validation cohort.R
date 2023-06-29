library(tidyverse)
library(DescTools)
library(ggplot2)
library(ggpubr)
library(reshape2)

Cancer<-readRDS("../ML_cancers_vs._cancerfree/Validation_Cancer_risk_score.rds")
Cancer<-Cancer[OrderMixed(names(Cancer))]
BRC<-readRDS("../ML_BRC_vs._others/Validation_BRC_risk_score.rds")
BRC<-BRC[OrderMixed(names(BRC))]
CRC<-readRDS("../ML_CRC_vs._others/Validation_CRC_risk_score.rds")
CRC<-CRC[OrderMixed(names(CRC))]
GC<-readRDS("../ML_GC_vs._others/Validation_GC_risk_score.rds")
GC<-GC[OrderMixed(names(GC))]
HCC<-readRDS("../ML_HCC_vs._others/Validation_HCC_risk_score.rds")
HCC<-HCC[OrderMixed(names(HCC))]
LC<-readRDS("../ML_LC_vs._others/Validation_LC_risk_score.rds")
LC<-LC[OrderMixed(names(LC))]

df<-data.frame(Cancer=Cancer,BRC=BRC,CRC=CRC,GC=GC,HCC=HCC,LC=LC)
df$actual=substr(rownames(df),1,3)
df$actual=gsub("_","",df$actual)
df$actual=gsub("NOR","Cancer-free individuals",df$actual)
df$actual=gsub("LC","Lung cancer patients",df$actual)
colnames(df)<-c("Cancer_207","BRC_39","CRC_70","GC_82","HCC_112","LC_51","actual")
write.csv(df,"validation cohort risk scores.csv")
df_plot<-melt(df)
df_plot$variable=gsub("Cancer_207","Cancer 207 cfRNA",df_plot$variable)
df_plot$variable=gsub("BRC_39","BRC 39 cfRNA",df_plot$variable)
df_plot$variable=gsub("CRC_70","CRC 70 cfRNA",df_plot$variable)
df_plot$variable=gsub("GC_82","GC 82 cfRNA",df_plot$variable)
df_plot$variable=gsub("HCC_112","HCC 112 cfRNA",df_plot$variable)
df_plot$variable=gsub("LC_51","LC 51 cfRNA",df_plot$variable)

df_plot$variable=factor(df_plot$variable,levels = c("LC 51 cfRNA","HCC 112 cfRNA","GC 82 cfRNA","CRC 70 cfRNA","BRC 39 cfRNA","Cancer 207 cfRNA"))

ggplot(df_plot, aes(x=variable, y=value,colour = actual)) + geom_boxplot(outlier.shape=NA)+ geom_point(position=position_jitterdodge())+
  scale_colour_manual(values=c("#4E62AB","#D6404E"))+labs(x="", y="Risk score")+ coord_flip()+  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic()+theme(text = element_text(size = 16), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.title = element_blank(),legend.position = "top")
ggsave(paste0("Validation cohort-median risk score.pdf"), width = 5, height = 5)
ggsave(paste0("Validation cohort-median risk score.png"), width = 5, height = 5, bg="white")
