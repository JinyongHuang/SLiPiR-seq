library(tidyverse)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(rstatix)
library(RColorBrewer)

meta_unique<-read.csv("../Metadata/Meta_with_clinical_and_Summary_all_samples.csv")

table(meta_unique$sample_type)

####Filter out bad samples
meta_good<-meta_unique[meta_unique$CleanRatio>20,]#bad library
meta_good<-meta_good[meta_good$Trimmed>2e6,]#low informative reads
meta_good<-meta_good[meta_good$rsRNA_ratio<30,]#rsRNA contemination
meta_good$gRNA_ratio<-meta_good$lncRNA_ratio+meta_good$mRNA_ratio
meta_good<-meta_good[meta_good$gRNA_ratio<30,]#cell contemination
meta_good<-meta_good[-(which(meta_good$sample_type=="LC_SZBU" & meta_good$AJCC.Stage=="IV")),]#one stage IV patient in the validation cohort
table(meta_good$sample_type)
write.csv(meta_good,"../Metadata/Meta_with_clinical_and_Summary_good_samples.csv",row.names = F)

####Used samples----
meta_use<-meta_good[meta_good$sample_type=="BRC_NXYK" | meta_good$sample_type=="CRC_NXYK" | meta_good$sample_type=="GC_NXYK" | meta_good$sample_type=="HCC_NXYK" |
                      meta_good$sample_type=="LC_SZDE" | meta_good$sample_type=="NOR_SZBA" | meta_good$sample_type=="LC_SZBU" | meta_good$sample_type=="NOR_SZDW",]
write.csv(meta_use,"../Metadata/Meta_with_clinical_and_Summary_used_samples.csv",row.names = F)

  
####box plot RNA ratio- lung cancer and healthy
meta_LCNOR<-meta_good[meta_good$sample_type=="LC_SZDE" | meta_good$sample_type=="NOR_SZBA" |
                      meta_good$sample_type=="LC_SZBU" | meta_good$sample_type=="NOR_SZDW",]
ratio_df<-meta_LCNOR[,c(3,grep("_ratio",colnames(meta_LCNOR)))]
ratio_df<-ratio_df[ , !(colnames(ratio_df) %in% "gRNA_ratio")]
ratio_df<-melt(ratio_df,"sample_type")
ratio_df$variable<-gsub("_ratio","",ratio_df$variable)
ratio_df$variable<-gsub("mRNA","msRNA",ratio_df$variable)
ratio_df$variable<-gsub("lncRNA","lsRNA",ratio_df$variable)
ratio_df$variable<-gsub("tRFs","tsRNA",ratio_df$variable)
ratio_df$sample_type<-factor(ratio_df$sample_type,levels = c("LC_SZDE","NOR_SZBA","LC_SZBU","NOR_SZDW"))
colnames(ratio_df)[2]<-"RNA" ##Dataframe shouldn't contain a column called "variable".
plot<-ggplot(ratio_df, aes(x=sample_type, y=value)) + geom_violin(aes(fill = sample_type))+ geom_boxplot(aes(fill = sample_type),width=0.2,color="black",outlier.shape=NA) + 
  facet_wrap(~RNA, scale="free")+ scale_fill_manual(values=c("#1F78B4","#A6CEE3","#1F78B4","#A6CEE3"))+ 
  labs(x="", y="Read counts (%)",title = "Read count ratio - all cases and controls")+
  theme_minimal()+ theme(legend.position="none",text = element_text(size=13),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot
ggsave("Read count ratio -LC cases and controls.pdf", width = 6, height = 8)
ggsave("Read count ratio -LC cases and controls.png", bg="white", width = 6, height = 8)

