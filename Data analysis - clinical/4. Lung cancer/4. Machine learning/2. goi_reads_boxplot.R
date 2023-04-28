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
goi_msRNA<-readLines("../LASSO/goi/20/msRNA name.txt")
goi_miRNA<-readLines("../Boruta/goi/20/miRNA name.txt")
goi_snRNA<-readLines("../LASSO/goi/20/snRNA name.txt")
goi_snoRNA<-readLines("../TopN/goi/snoRNA Top N selected name.txt")
goi_tsRNA<-readLines("../LASSO/goi/20/tsRNA name.txt")
goi_list<-list(msRNA=goi_msRNA,miRNA=goi_miRNA,
               tsRNA=goi_tsRNA,snRNA=goi_snRNA,snoRNA=goi_snoRNA,
               cfRNA=c(goi_msRNA,goi_miRNA,goi_tsRNA,goi_snoRNA))

####case and control----
##each RNA
subset<-AllRNA[GOI_lung,]
GOI_lung<-c(goi_msRNA,goi_miRNA,goi_tsRNA,goi_snRNA,goi_snoRNA)
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

combine<-ggarrange(plotlist[["AP1B1"]],plotlist[["FOXF1"]],plotlist[["TMC5"]],
                   plotlist[["hsa-miR-25-3p"]],plotlist[["hsa-miR-1180-3p"]],plotlist[["hsa-miR-3960"]],
                   plotlist[["SNORD6"]],plotlist[["SNORD25"]],plotlist[["SNORD53"]],
                   plotlist[["tRF-36-PS5P4PW3FJHPEZE"]],plotlist[["tRF-39-18VBY9P9KHM2ORIZ"]],plotlist[["tRF-45-Z3R918VBY9PYKHM26R"]],
                   ncol = 3, nrow = 4)
ggsave("Combined 12.pdf", width=10.5, height=14)
ggsave("Combined 12.png", width=10.5, height=14, units="in", bg="white")


##cumulative RNA
dir.create(paste0(getwd(),"/Boxplot_cumulative case and control"), showWarnings = FALSE)
LCNOR<-meta[(meta$sample_type=="LC_SZDE" | meta$sample_type=="NOR_SZBA" | meta$sample_type=="LC_SZBU" | meta$sample_type=="NOR_SZDW"),]
for (i in names(goi_list)) {
  subset<-AllRNA[goi_list[[i]],LCNOR$ID]
  subset<-data.frame(t(subset),check.names = FALSE)
  subset<-as.data.frame(apply(subset, 1, sum))
  colnames(subset)<-"RNA"
  subset$Group<-LCNOR$Group
  subset$Type<-LCNOR$sample_type
  subset$Type<-gsub("LC_SZDE","LC-Disc",subset$Type)
  subset$Type<-gsub("NOR_SZBA","NOR-Disc",subset$Type)
  subset$Type<-gsub("LC_SZBU","LC-Vali",subset$Type)
  subset$Type<-gsub("NOR_SZDW","NOR-Vali",subset$Type)
  subset$Type<-factor(subset$Type,levels = c("LC-Disc","NOR-Disc","LC-Vali","NOR-Vali"))
  stat.test <- subset %>%
    t_test(RNA ~ Type, paired = FALSE, var.equal = FALSE, p.adjust.method = "none")%>%
    add_significance("p")%>%
    add_xy_position(fun="max")
  stat.test<-stat.test[c(1,6),]
  stat.test$y.position<-mean(stat.test$y.position)
  if (i=="cfRNA") {
    plot <- ggplot(subset, aes(x=Type, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+ scale_fill_manual(values=c("#4E62AB","#D6404E"))+ 
      geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.8)+
      labs(title=paste0(i),x="", y = expression('Cumulative log'[2]*'(RPM)'))+ theme_minimal()+
      theme(axis.title = element_text(size = 14), axis.text = element_text(colour = "black",size = 12),
            plot.title = element_text(hjust=0.5,face="bold"),legend.position = "none")}
  if (i!="cfRNA") {
    plot <- ggplot(subset, aes(x=Type, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+ scale_fill_manual(values=c("#4E62AB","#D6404E"))+ 
      geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.8)+
      labs(title=paste0(i), x="", y = expression('Cumulative log'[2]*'(RPM)'))+ theme_minimal()+
      theme(axis.title = element_text(size = 14), axis.text = element_text(colour = "black",size = 12), 
            plot.title = element_text(hjust=0.5,face="bold"),legend.position = "none")}
  plot + stat_pvalue_manual(stat.test, label = "p.signif",tip.length = 0.02,bracket.shorten = 0.1)
  ggsave(paste0("Boxplot_cumulative case and control/",i,".png"),device = "png",bg="white",width = 4, height = 4)
  ggsave(paste0("Boxplot_cumulative case and control/",i,".pdf"),device = "png",width = 4, height = 4)
}


