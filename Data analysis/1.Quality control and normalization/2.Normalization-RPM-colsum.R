library(tidyverse)
library(DescTools)
#library(edgeR) #cpm()
library(gtools)
library(ggplot2)
library(ggpubr)
####meta----
meta<-read.csv("../Metadata/Meta_with_clinical_and_Summary_used_samples_20230829.csv", header = T, row.names = 1)
table(meta$sample_type)
meta<-meta[OrderMixed(meta$ID),]
####RNA matrix----
###mRNA
mRNA<-read.csv("../Raw read counts/mRNA.csv", header = T, row.names = 1)
mRNA<-mRNA[,meta$ID]
###lncRNA
lncRNA<-read.csv("../Raw read counts/lncRNA.csv", header = T, row.names = 1)
lncRNA<-lncRNA[,meta$ID]
###miRNA
miRNA<-read.csv("../Raw read counts/miRNA.csv", header = T, row.names = 1)
miRNA<-miRNA[,meta$ID]
###piRNA
piRNA<-read.csv("../Raw read counts/piRNA.csv", header = T, row.names = 1)
piRNA<-piRNA[,meta$ID]
###tRFs
tsRNA<-read.csv("../Raw read counts/tRFs.csv", header = T, row.names = 1)
tsRNA<-tsRNA[,meta$ID]
###rsRNA
rsRNA<-read.csv("../Raw read counts/rsRNA.csv", header = T, row.names = 1)
rsRNA<-rsRNA[,meta$ID]
###ysRNA
ysRNA<-read.csv("../Raw read counts/ysRNA.csv", header = T, row.names = 1)
ysRNA<-ysRNA[,meta$ID]
###snRNA
snRNA<-read.csv("../Raw read counts/snRNA.csv", header = T, row.names = 1)
snRNA<-snRNA[,meta$ID]
###snoRNA
snoRNA<-read.csv("../Raw read counts/snoRNA.csv", header = T, row.names = 1)
snoRNA<-snoRNA[,meta$ID]
##All RNA
cfRNA<-rbind(mRNA,lncRNA,miRNA,piRNA,snRNA,snoRNA,tsRNA,rsRNA,ysRNA)
cfRNA<-cfRNA[,OrderMixed(colnames(cfRNA))]
colsum=colSums(cfRNA)
###List
RNA<-list(mRNA=mRNA, lncRNA=lncRNA, miRNA=miRNA, piRNA=piRNA, snRNA=snRNA,snoRNA=snoRNA, tsRNA=tsRNA, rsRNA=rsRNA, ysRNA=ysRNA,cfRNA=cfRNA)

####log2(cpm+1) normalization----
##Normalize by studied RNA reads
RNA_cpm<-list()
RNA_log2cpm<-list()
for (i in 1:length(RNA)) {
  print(paste0("Start working on ",names(RNA)[i]))
  subset <- RNA[[i]]
  subset<-subset[,OrderMixed(colnames(subset))]
  print(table(colnames(subset)==meta$ID))
  print(table(colnames(subset)==names(colsum)))
  CPM<-apply(subset,1,function(x) x/colsum*1e+6)
  CPM<-as.data.frame(t(CPM))
  write.table(CPM, paste0("../Normalization_RPM_colsum/",names(RNA)[i],"-cpm.txt"),sep = "\t", col.names = NA, quote = FALSE)
  RNA_cpm[[i]]<-CPM
  CPM[]<-sapply(CPM, function(x) log2(x+1))
  write.table(CPM, paste0("../Normalization_RPM_colsum/",names(RNA)[i],"-log2(cpm).txt"),sep = "\t", col.names = NA, quote = FALSE)
  RNA_log2cpm[[i]]<-CPM
}
names(RNA_cpm)=names(RNA)
names(RNA_log2cpm)=names(RNA)

####Pearson's correlation coefficient (R) & Coefficient of Determination (R-Squared)----
meta<-meta[meta$sample_type=="BRC_NXYK" | meta$sample_type=="CRC_NXYK" | meta$sample_type=="GC_NXYK" | meta$sample_type=="HCC_NXYK" |
             meta$sample_type=="LC_SZDE" | meta$sample_type=="NOR_SZBA",]
combination<-data.frame(t(combn(unique(meta$sample_type),2)))
combination$X3=paste(combination$X1,"vs.",combination$X2)
corcoef <-matrix(nrow=nrow(combination),ncol=length(RNA))
rownames(corcoef)<-str_c(combination[,1]," vs. ",combination[,2])
colnames(corcoef)<-names(RNA)
RSquared<-matrix(nrow=nrow(combination),ncol=length(RNA))
rownames(RSquared)<-str_c(combination[,1]," vs. ",combination[,2])
colnames(RSquared)<-names(RNA)
for (i in 1:nrow(combination)) {
  print(combination[i,3])
  group1<-combination[i,1]
  id1<-meta[meta$sample_type==group1,][,"ID"]
  group2<-combination[i,2]
  id2<-meta[meta$sample_type==group2,][,"ID"]
  for (r in 1:length(RNA)){
    temp1<-RNA_log2cpm[[r]][,id1]
    temp2<-RNA_log2cpm[[r]][,id2]
    temp <- merge(apply(temp1,1,mean),apply(temp2,1,mean),by=0)
    corcoef[i,r]=cor(temp$x,temp$y)
    model <- lm(y~x, data=temp)
    RSquared[i,r]=summary(model)$r.squared
  }
}
write.csv(corcoef,"Correlation coefficient_by cohort.csv")
write.csv(RSquared,"Coefficient of Determination (R-Squared) by cohort.csv")

corcoef_df<-data.frame(corcoef)
rownames(corcoef_df)<-gsub("BRC_NXYK","BRC",rownames(corcoef_df))
rownames(corcoef_df)<-gsub("CRC_NXYK","CRC",rownames(corcoef_df))
rownames(corcoef_df)<-gsub("GC_NXYK","GC",rownames(corcoef_df))
rownames(corcoef_df)<-gsub("HCC_NXYK","HCC",rownames(corcoef_df))
rownames(corcoef_df)<-gsub("LC_SZDE","LC",rownames(corcoef_df))
rownames(corcoef_df)<-gsub("NOR_SZBA","NOR",rownames(corcoef_df))
write.csv(corcoef_df,"Correlation coefficient_by cohort renamed.csv")

####Scatterplot for publication----
plot_list<-list()
group1<-"LC_SZDE"
id1<-meta[meta$sample_type==group1,][,"ID"]
group2<-"NOR_SZBA"
id2<-meta[meta$sample_type==group2,][,"ID"]

for (r in names(RNA)[1:9]){
  temp1<-RNA_log2cpm[[r]][,id1]
  temp2<-RNA_log2cpm[[r]][,id2]
  temp <- merge(apply(temp1,1,mean),apply(temp2,1,mean),by=0)
  plot <- ggplot(temp, aes(x=x, y=y)) + geom_point(size=1,alpha=0.8,color="#D6404E") + geom_rug() + theme_classic() +
    stat_cor(aes(label = paste(after_stat(r.label))),digits = 3,method = "pearson", label.x = range(temp1[,1])[2]*0.6, label.y = range(temp1[,1])[2]*0.15, size=4)+
    labs(title=paste0(r)) +theme(axis.title=element_blank(),axis.text = element_text(color = "black"),
                                 plot.title = element_text(hjust=0.5,face = "bold"))
  
  plot_list[[r]]<-plot
}
plot<-ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
          plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],ncol = 3, nrow = 3)
annotate_figure(plot,
                left = text_grob(expression("Log"[2]*"(RPM) of controls"),rot = 90),
                bottom=text_grob(expression("Log"[2]*"(RPM) of cases")))
ggsave("Pearson's correlation 9 RNA.pdf", width = 9, height = 9)
ggsave("Pearson's correlation 9 RNA.png", width = 9, height = 9, units = "in", bg="white")

temp1<-RNA_log2cpm[["cfRNA"]][,id1]
temp2<-RNA_log2cpm[["cfRNA"]][,id2]
temp <- merge(apply(temp1,1,mean),apply(temp2,1,mean),by=0)
plot <- ggplot(temp, aes(x=x, y=y)) + geom_point(size=0.3,alpha=0.8,color="#D6404E") + geom_rug() + theme_minimal() +
  stat_cor(aes(label = after_stat(r.label)),digits = 3,method = "pearson", label.x = 8, label.y = 2.5, size=2)+
  stat_cor(aes(label = after_stat(p.label)),digits = 3,method = "pearson", label.x = 8, label.y = 1, size=2)+
  labs(x=expression("Log"[2]*"(RPM) of cases"), y = expression("Log"[2]*"(RPM) of controls")) + theme_classic()+
  theme(axis.text = element_text(color = "black",size=6),axis.title = element_text(size=6),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),axis.line = element_line(linewidth=0.15))
ggsave("Pearson's correlation all cfRNA.pdf", width = 5, height = 5, units="cm")
ggsave("Pearson's correlation all cfRNA.tiff", width = 5, height = 5, units = "cm", bg="white")

