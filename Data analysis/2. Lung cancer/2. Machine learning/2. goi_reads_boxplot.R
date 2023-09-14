library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(reshape2)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
####Sample ID----
meta<-read.csv("../../../Metadata/Meta_with_clinical_and_Summary_good_samples.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="NOR_SZBA" | sample_type=="LC_SZBU" | sample_type=="NOR_SZDW" |
                 sample_type=="BRC_NXYK" | sample_type=="CRC_NXYK" | sample_type=="GC_NXYK" | sample_type=="HCC_NXYK")

####RNA cpm
AllRNA<-read.table("../../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[,meta$ID]

##GOI
goi_mRNA<-readLines("../LASSO_10/goi/mRNA name.txt")
goi_miRNA<-readLines("../LASSO_10/goi/miRNA name.txt")
goi_snRNA<-readLines("../LASSO_10/goi/snRNA name.txt")
goi_snoRNA<-readLines("../LASSO_10/goi/snoRNA name.txt")
goi_tsRNA<-readLines("../LASSO_10/goi/tsRNA name.txt")
goi_rsRNA<-readLines("../LASSO_10/goi/rsRNA name.txt")
goi_ysRNA<-readLines("../LASSO_10/goi/ysRNA name.txt")
goi_list<-list(mRNA=goi_mRNA,miRNA=goi_miRNA,
               tsRNA=goi_tsRNA,snRNA=goi_snRNA,snoRNA=goi_snoRNA,
               rsRNA=goi_rsRNA,ysRNA=goi_ysRNA)

####case and control----
##each RNAsubset<-AllRNA[GOI_lung,]
GOI_lung<-c(goi_mRNA,goi_miRNA,goi_tsRNA,goi_snRNA,goi_snoRNA,goi_rsRNA,goi_ysRNA)
dir.create(paste0(getwd(),"/Boxplot_each case and control"), showWarnings = FALSE)
subset<-AllRNA[GOI_lung,]
subset<-data.frame(t(subset),check.names = FALSE)
subset$Group<-matrix(unlist(str_split(rownames(subset),"_")),ncol=nrow(subset),nrow=3)[1,]
subset<-subset[subset$Group=="LC" | subset$Group=="NOR",]
subset$Type<-matrix(unlist(str_split(rownames(subset),"_")),ncol=nrow(subset),nrow=3)[2,]
subset$Type<-gsub("SZDE","LC-Disc",subset$Type)
subset$Type<-gsub("SZBA","NOR-Disc",subset$Type)
subset$Type<-gsub("SZBU","LC-Vali",subset$Type)
subset$Type<-gsub("SZDW","NOR-Vali",subset$Type)
subset$Type<-factor(subset$Type,levels = c("LC-Disc","NOR-Disc","LC-Vali","NOR-Vali"))
N=1:(ncol(subset)-2)
plotlist<-list()
for (i in N) {
  df<-subset[,c(i,ncol(subset)-1,ncol(subset))]
  colnames(df)[1]<-"RNA"
  stat.test <- df %>%
    t_test(RNA ~ Type, paired = FALSE, var.equal = FALSE, p.adjust.method = "none")%>%
    add_significance("p")%>%
    add_xy_position(fun="max")
  stat.test<-stat.test[c(1,6),]
  stat.test$y.position<-mean(stat.test$y.position)
  if (stat.test$p[1]<0.05 & stat.test$p[2]<0.05) {
    plot <- ggplot(df, aes(x=Type, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+ scale_fill_manual(values=c("#D6404E","#4E62AB"))+ 
    geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.8)+
    labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*'(RPM)'))+ theme_minimal()+
    theme(axis.title = element_text(size = 14), axis.text = element_text(colour = "black",size = 12), 
          plot.title = element_text(hjust=0.5,face="bold"),legend.position = "none")
    plot <- plot + stat_pvalue_manual(stat.test, label = "p.signif", bracket.shorten = 0.5)
    plotlist[[colnames(subset)[i]]]<-plot
  ggsave(paste0("Boxplot_each case and control/",colnames(subset)[i],".png"),device = "png",bg="white",width = 4, height = 4)
  ggsave(paste0("Boxplot_each case and control/",colnames(subset)[i],".pdf"),device = "pdf",width = 4, height = 4)
}}

combine<-ggarrange(plotlist[["AP1B1"]],plotlist[["FOXF1"]],plotlist[["RPL23A"]],plotlist[["TMC5"]],
                   plotlist[["hsa-miR-185-5p"]],plotlist[["hsa-miR-1180-3p"]],plotlist[["hsa-miR-3960"]],
                   plotlist[["SNORD16"]],plotlist[["SNORD22"]],plotlist[["SNORD53"]],
                   plotlist[["tRF-38-HMI8W47W1R7HFEV"]],plotlist[["tRF-33-R008R959KUMK9"]],plotlist[["tRF-23-BZBZOS4YV"]],
                   plotlist[["RNU4-46P"]],plotlist[["Homo-28S-2571"]],
                   ncol = 3, nrow = 5)
ggsave("Combined 15.pdf", width=10.5, height=18)
ggsave("Combined 15.png", width=10.5, height=18, units="in", bg="white")
#ggsave("Combined 15.tiff", width=10.5, height=18, units="in", bg="white")

