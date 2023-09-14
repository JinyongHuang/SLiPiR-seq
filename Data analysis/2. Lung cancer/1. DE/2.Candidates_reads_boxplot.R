library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(reshape2)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(RColorBrewer)

####Sample ID----
meta<-read.csv("../../Metadata/Meta_with_clinical_and_Summary_good_samples.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="NOR_SZBA" | sample_type=="LC_SZBU" | sample_type=="NOR_SZDW")
table(meta$sample_type)

####RNA cpm
AllRNA<-read.table("../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[,meta$ID]

####Significant RNA
RNA_name <- list.files(path = "../../Normalization_RPM_colsum/", pattern = "cpm.txt$", full.names = TRUE, recursive = TRUE) 
RNA_name <- RNA_name[grep("cfRNA",RNA_name, invert = TRUE)]
RNA_name <- lapply(RNA_name, read.table, header=TRUE, sep="\t",row.names=1)
RNA_name <- lapply(RNA_name, row.names)
names(RNA_name)<- c("lncRNA",  "miRNA",  "mRNA",  "piRNA",  "rsRNA",  "snoRNA", "snRNA",  "tsRNA",  "ysRNA")
sig_cfRNA<-read.table("../DE_LC_SZDE_vs_NOR_SZBA/DESeq2 candidate results-cfRNA-Lung cancer vs. Healthy.txt", header = T)
sig_cfRNA<-sig_cfRNA[sig_cfRNA$log2FoldChange>0,]
sig_list<-list()
Distribution<-double()
for (RNA in names(RNA_name)) {
  sig_list[[RNA]]<-rownames(sig_cfRNA)[(rownames(sig_cfRNA) %in% RNA_name[[RNA]])]
  Distribution[RNA]=table(rownames(sig_cfRNA) %in% RNA_name[[RNA]])["TRUE"]
}
sum(Distribution)

####early and late stage----
discovery<-meta[(meta$sample_type=="LC_SZDE" | meta$sample_type=="NOR_SZBA"),]
##each RNA
dir.create(paste0(getwd(),"/Boxplot_each early and late stage"), showWarnings = FALSE)
subset<-AllRNA[rownames(sig_cfRNA),discovery$ID]
subset<-data.frame(t(subset),check.names = FALSE)
subset$Group<-discovery$Group
subset$Type<-discovery$Stage
subset$Type[is.na(subset$Type)]<-"Cancer-free"
subset$Type<-gsub("Early","Early stage",subset$Type)
subset$Type<-gsub("Late","Late stage",subset$Type)
subset$Type<-factor(subset$Type,levels = c("Cancer-free","Early stage","Late stage"))
N=1:(ncol(subset)-2)
registerDoParallel(cl<-makeCluster(10))
result_seed<-foreach(i=N, .combine="c",.packages=c("tidyverse","ggplot2","ggpubr","rstatix","reshape2")) %dopar% {
  df<-subset[,c(i,ncol(subset)-1,ncol(subset))]
  colnames(df)[1]<-"RNA"
  stat.test <- df %>%
    t_test(RNA ~ Type, paired = FALSE, var.equal = FALSE, p.adjust.method = "none", ref.group = "Cancer-free")%>%
    add_significance("p")%>%
    add_xy_position(fun="max")
  if (stat.test$p[1]<0.05 & stat.test$p[2]<0.05) {
    plot <- ggplot(df, aes(x=Type, y=RNA)) + geom_boxplot(aes(fill = Type),outlier.shape=NA)+ scale_fill_manual(values=c("#4E62AB","#D6404E","#9E0142"))+ 
      geom_jitter(shape=16, position=position_jitter(0.2))+
      labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*' (RPM+1)'))+ theme_minimal()+
      theme(text = element_text(size = 14),legend.position = "none")
    plot + stat_pvalue_manual(stat.test, label = "p.signif",tip.length = 0)
    ggsave(paste0("Boxplot_each early and late stage/",colnames(subset)[i],".png"),device = "png",bg="white",width = 6, height = 6)
  }}
stopCluster(cl)
  
##cumulative RNA
dir.create(paste0(getwd(),"/Boxplot_cumulative early and late stage"), showWarnings = FALSE)
plot_list<-list()
for (RNA in names(sig_list)) {
  subset<-AllRNA[sig_list[[RNA]],discovery$ID]
  subset<-data.frame(t(subset),check.names = FALSE)
  subset<-as.data.frame(apply(subset, 1, sum))
  colnames(subset)<-"RNA"
  subset$Group<-discovery$Group
  subset$Type<-discovery$Stage
  subset$Type[is.na(subset$Type)]<-"Cancer-free"
  subset$Type<-gsub("Early","Early stage",subset$Type)
  subset$Type<-gsub("Late","Late stage",subset$Type)
  subset$Type<-factor(subset$Type,levels = c("Cancer-free","Early stage","Late stage"))
  stat.test <- subset %>%
    t_test(RNA ~ Type, paired = FALSE, var.equal = FALSE, p.adjust.method = "none")%>%
    add_significance("p")%>%
    add_xy_position(fun="max")
  max<-max(subset$RNA)
  stat.test$y.position[1]<-max+(max/100*2)
  stat.test$y.position[2]<-max+(max/100*7)
  stat.test$y.position[3]<-max+(max/100*2)
  plot <- ggplot(subset, aes(x=Type, y=RNA)) + geom_boxplot(aes(fill = Type),outlier.shape=NA)+ scale_fill_manual(values=c("#4E62AB","#D6404E","#9E0142"))+ 
    geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.8)+
    labs(title=paste0(RNA))+ theme_classic()+
    theme(text = element_text(size = 13),legend.position = "none",plot.title = element_text(hjust=0.5,face="bold"),
          axis.text = element_text(colour = "black"),axis.title = element_blank())
  plot <- plot + stat_pvalue_manual(stat.test, label = "p.signif",tip.length = 0.02,bracket.shorten = 0.1)
  ggsave(paste0("Boxplot_cumulative early and late stage/",RNA,".png"),device = "png",bg="white",width = 4, height = 4)
  plot_list[[RNA]]<-plot
}

combine<-ggarrange(plot_list[[3]], plot_list[[1]], plot_list[[2]],
          plot_list[[4]], plot_list[[7]], plot_list[[6]], 
          plot_list[[8]], plot_list[[5]], plot_list[[9]], 
          ncol = 3, nrow = 3)
annotate_figure(combine,
                left = text_grob(expression('Cumulative log'[2]*'(RPM)'),rot = 90,size=16))
ggsave("Boxplot_cumulative - combined 9.pdf", width=9, height=9)
ggsave("Boxplot_cumulative - combined 9.png", width=9, height=9, units="in", bg="white")

##cfRNA
subset<-AllRNA[rownames(sig_cfRNA),discovery$ID]
subset<-data.frame(t(subset),check.names = FALSE)
subset<-as.data.frame(apply(subset, 1, sum))
colnames(subset)<-"RNA"
subset$Group<-discovery$Group
subset$Type<-discovery$Stage
subset$Type[is.na(subset$Type)]<-"Cancer-free"
subset$Type<-gsub("Early","Early stage",subset$Type)
subset$Type<-gsub("Late","Late stage",subset$Type)
subset$Type<-factor(subset$Type,levels = c("Cancer-free","Early stage","Late stage"))
stat.test <- subset %>%
  t_test(RNA ~ Type, paired = FALSE, var.equal = FALSE, p.adjust.method = "none")%>%
  add_significance("p")%>%
  add_xy_position(fun="max")
max<-max(subset$RNA)
stat.test$y.position[1]<-max+(max/100*2)
stat.test$y.position[2]<-max+(max/100*7)
stat.test$y.position[3]<-max+(max/100*2)
plot <- ggplot(subset, aes(x=Type, y=RNA)) + geom_boxplot(aes(fill = Type),outlier.shape=NA)+ scale_fill_manual(values=c("#4E62AB","#D6404E","#9E0142"))+ 
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.8, size=0.5)+
  labs(x="", y = expression('Cumulative log'[2]*'(RPM)'))+ theme_classic()+
  theme(text = element_text(size = 6),legend.position = "none", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),axis.line = element_line(linewidth = 0.15),
        axis.text = element_text(colour = "black"),axis.text.x=element_text(angle = 30,hjust=1))
plot + stat_pvalue_manual(stat.test, label = "p.signif",tip.length = 0.02,bracket.shorten = 0.1, size=2)
ggsave(paste0("Boxplot_cumulative early and late stage cfRNA.png"),device = "png",bg="white",width = 4, height = 5, units="cm")
ggsave(paste0("Boxplot_cumulative early and late stage cfRNA.pdf"),device = "pdf",width = 4, height = 5, units="cm")
