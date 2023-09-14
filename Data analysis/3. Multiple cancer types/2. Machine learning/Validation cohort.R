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
colnames(df)<-c("Cancer_65","BRC_21","CRC_33","GC_36","HCC_33","LC_30","actual")
write.csv(df,"validation cohort risk scores.csv")
df_plot<-melt(df)
df_plot$variable=gsub("Cancer_65","Cancer 65 cfRNAs",df_plot$variable)
df_plot$variable=gsub("BRC_21","BRC 21 cfRNAs",df_plot$variable)
df_plot$variable=gsub("CRC_33","CRC 33 cfRNAs",df_plot$variable)
df_plot$variable=gsub("GC_36","GC 36 cfRNAs",df_plot$variable)
df_plot$variable=gsub("HCC_33","HCC 33 cfRNAs",df_plot$variable)
df_plot$variable=gsub("LC_30","LC 30 cfRNAs",df_plot$variable)


df_plot$variable=factor(df_plot$variable,levels = c("LC 30 cfRNAs","HCC 33 cfRNAs","GC 36 cfRNAs","CRC 33 cfRNAs","BRC 21 cfRNAs","Cancer 65 cfRNAs"))

ggplot(df_plot, aes(x=variable, y=value,colour = actual)) + geom_boxplot(outlier.shape=NA)+ geom_point(position=position_jitterdodge())+
  scale_colour_manual(values=c("#4E62AB","#D6404E"))+labs(x="", y="Risk score")+ coord_flip()+  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  theme_classic()+theme(text = element_text(size = 16), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),legend.title = element_blank(),legend.position = "top")
ggsave(paste0("Validation cohort-median risk score.pdf"), width = 5, height = 5)
ggsave(paste0("Validation cohort-median risk score.png"), width = 5, height = 5, bg="white")


####Discovery cohort risk scores
Disc_Cancer<-readRDS("../ML_cancers_vs._cancerfree/Cancers-Discovery cohort risk scores.rds")
Disc_Cancer<-Disc_Cancer[OrderMixed(rownames(Disc_Cancer)),]
colnames(Disc_Cancer)[1]<-'Cancer panel'
Disc_BRC<-readRDS("../ML_BRC_vs._others/BRC-Discovery cohort risk scores.rds")
Disc_BRC<-Disc_BRC[OrderMixed(rownames(Disc_BRC)),]
colnames(Disc_BRC)[1]<-'BRC panel'
Disc_BRC<-Disc_BRC[1]
Disc_CRC<-readRDS("../ML_CRC_vs._others/CRC-Discovery cohort risk scores.rds")
Disc_CRC<-Disc_CRC[OrderMixed(rownames(Disc_CRC)),]
colnames(Disc_CRC)[1]<-'CRC panel'
Disc_CRC<-Disc_CRC[1]
Disc_GC<-readRDS("../ML_GC_vs._others/GC-Discovery cohort risk scores.rds")
Disc_GC<-Disc_GC[OrderMixed(rownames(Disc_GC)),]
colnames(Disc_GC)[1]<-'GC panel'
Disc_GC<-Disc_GC[1]
Disc_HCC<-readRDS("../ML_HCC_vs._others/HCC-Discovery cohort risk scores.rds")
Disc_HCC<-Disc_HCC[OrderMixed(rownames(Disc_HCC)),]
colnames(Disc_HCC)[1]<-'HCC panel'
Disc_HCC<-Disc_HCC[1]
Disc_LC<-readRDS("../ML_LC_vs._others/LC-Discovery cohort risk scores.rds")
Disc_LC<-Disc_LC[OrderMixed(rownames(Disc_LC)),]
colnames(Disc_LC)[1]<-'LC panel'
Disc_LC<-Disc_LC[1]

Disc_combine<-cbind(Disc_Cancer,Disc_BRC,Disc_CRC,Disc_GC,Disc_HCC,Disc_LC)
write.csv(Disc_combine,"Discovery cohort risk scores.csv")
