library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggrepel)
library(RColorBrewer)
library(patchwork)#Combine ggplot
library(DESeq2)#Differential analysis
library(EnhancedVolcano)#Volcano
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

####Sample ID----
meta<-read.csv("../../Metadata/Meta_with_clinical_and_Summary_good_samples.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="NOR_SZBA" |
                 sample_type=="BRC_NXYK" | sample_type=="CRC_NXYK" | sample_type=="GC_NXYK" | sample_type=="HCC_NXYK")
meta$group<-ifelse(meta$Group=="Liver cancer","Cases","Controls")
meta$group<-factor(meta$group)
table(meta$sample_type)
case="Cases"
control="Controls"
ncase=table(meta$group==case)["TRUE"]
ncontrol=table(meta$group==control)["TRUE"]

####RNA matrix----
##mRNA
msRNA<-read.csv("../../Raw read counts/mRNA.csv", header = T, row.names = 1)
msRNA<-msRNA[,meta$ID]
##miRNA
miRNA<-read.csv("../../Raw read counts/miRNA.csv", header = T, row.names = 1)
miRNA<-miRNA[,meta$ID]
##lncRNA
lsRNA<-read.csv("../../Raw read counts/lncRNA.csv", header = T, row.names = 1)
lsRNA<-lsRNA[,meta$ID]
##piRNA
piRNA<-read.csv("../../Raw read counts/piRNA.csv", header = T, row.names = 1)
piRNA<-piRNA[,meta$ID]
##tRFs
tsRNA<-read.csv("../../Raw read counts/tRFS.csv", header = T, row.names = 1)
tsRNA<-tsRNA[,meta$ID]
##rsRNA
rsRNA<-read.csv("../../Raw read counts/rsRNA.csv", header = T, row.names = 1)
rsRNA<-rsRNA[,meta$ID]
##ysRNA
ysRNA<-read.csv("../../Raw read counts/ysRNA.csv", header = T, row.names = 1)
ysRNA<-ysRNA[,meta$ID]
##snRNA
snRNA<-read.csv("../../Raw read counts/snRNA.csv", header = T, row.names = 1)
snRNA<-snRNA[,meta$ID]
##snoRNA
snoRNA<-read.csv("../../Raw read counts/snoRNA.csv", header = T, row.names = 1)
snoRNA<-snoRNA[,meta$ID]
##All RNA
cfRNA<-rbind(msRNA,lsRNA,miRNA,piRNA,tsRNA,rsRNA,ysRNA,snRNA,snoRNA)
RNA_name<-list(msRNA=msRNA,lsRNA=lsRNA,miRNA=miRNA,piRNA=piRNA,tsRNA=tsRNA,rsRNA=rsRNA,ysRNA=ysRNA,snRNA=snRNA,snoRNA=snoRNA)
RNA_name <- lapply(RNA_name, row.names)

##DESeq2----
dds<-DESeqDataSetFromMatrix(cfRNA, meta, ~ group)
dds<-DESeq(dds)   
results<-results(dds, contrast = c("group", case, control)) 
summary(results)
results<-data.frame(results)
results<-results[!(is.na(results$padj)),]
results<-results[order(results$pvalue),]
write.table(results, paste0("DESeq2 results-cfRNA-",case," vs. ",control,".txt"),sep = "\t", col.names = NA, quote = FALSE)
sig_results<-results[results$padj < 0.1,]
write.table(sig_results, paste0("DESeq2 significant results-cfRNA-",case," vs. ",control,".txt"),sep = "\t", col.names = NA, quote = FALSE)
candi_results<-results[results$padj < 0.1 & results$baseMean>=10 & results$log2FoldChange>=1,]
write.table(candi_results, paste0("DESeq2 candidate results-cfRNA-",case," vs. ",control,".txt"),sep = "\t", col.names = NA, quote = FALSE)

##Distribution
all_list<-list()
all_name<-list()
all_distribution<-double()
for (RNA in names(RNA_name)) {
  all_list[[RNA]]<-results[(rownames(results) %in% RNA_name[[RNA]]),]
  all_name[[RNA]]<-rownames(results)[(rownames(results) %in% RNA_name[[RNA]])]
  all_distribution[RNA]=table(rownames(results) %in% RNA_name[[RNA]])["TRUE"]
}
sig_list<-list()
sig_name<-list()
sig_distribution<-double()
for (RNA in names(RNA_name)) {
  sig_list[[RNA]]<-sig_results[(rownames(sig_results) %in% RNA_name[[RNA]]),]
  sig_name[[RNA]]<-rownames(sig_results)[(rownames(sig_results) %in% RNA_name[[RNA]])]
  sig_distribution[RNA]=table(rownames(sig_results) %in% RNA_name[[RNA]])["TRUE"]
}
candi_list<-list()
candi_name<-list()
candi_distribution<-double()
for (RNA in names(RNA_name)) {
  candi_list[[RNA]]<-candi_results[(rownames(candi_results) %in% RNA_name[[RNA]]),]
  candi_name[[RNA]]<-rownames(candi_results)[(rownames(candi_results) %in% RNA_name[[RNA]])]
  candi_distribution[RNA]=table(rownames(candi_results) %in% RNA_name[[RNA]])["TRUE"]
}

####Summary
DESeq_result<-matrix(ncol=4,nrow=9)
colnames(DESeq_result)<-c("gain","loss","gain_filter","loss_filter")
rownames(DESeq_result)<-names(RNA_name)
for (RNA in names(RNA_name)) {
  temp<-all_list[[RNA]]
  gain<-nrow(subset(temp, padj<0.1 & log2FoldChange>0))
  DESeq_result[RNA,"gain"]<-gain
  loss<-nrow(subset(temp, padj<0.1 & log2FoldChange<0))
  DESeq_result[RNA,"loss"]<-loss
  gain_filter<-nrow(temp[temp$padj < 0.1 & temp$baseMean>=10 & temp$log2FoldChange>=1,])
  DESeq_result[RNA,"gain_filter"]<-gain_filter
  loss_filter<-nrow(temp[temp$padj < 0.1 & temp$baseMean>=10 & temp$log2FoldChange<=(-1),])
  DESeq_result[RNA,"loss_filter"]<-loss_filter
}
DESeq_result<-data.frame(DESeq_result)
sum<-apply(DESeq_result, 2, sum)
DESeq_result<-rbind(DESeq_result,sum)
rownames(DESeq_result)[10]<-"cfRNA"
write.csv(DESeq_result,"gain and loss.csv")

##Volcano plot----
plot<-list()
for (RNA in names(RNA_name)) {
  temp<-all_list[[RNA]]
  x1=min(temp$log2FoldChange)
  if (x1>(-2)) {x1=(-2)}
  if (x1<=(-2)) {x1=x1}
  x2=max(temp$log2FoldChange)
  y=-log10(min(temp$padj))
  increase<-nrow(subset(temp, padj<0.1 & log2FoldChange>0))
  decrease<-nrow(subset(temp, padj<0.1 & log2FoldChange<0))
  temp_Volcano<-EnhancedVolcano(temp, lab=NA,
                                x = "log2FoldChange", y = "padj", 
                                title = paste0(RNA), subtitle = paste0("loss ",decrease,"  |  gain ",increase),
                                pCutoff = 0.1, FCcutoff = 1, xlim = c(x1,x2), ylim = c(0,y), pointSize = 1, legendPosition = 'none',
                                col = c("black","black","#4E62AB","#D6404E"), colAlpha = 0.9, legendLabels = c("NS", "FC", "FDR", "FDR & FC"))
  temp_Volcano<-temp_Volcano + theme(axis.title = element_blank(),axis.text = element_text(color = "black"),
                                     plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
  plot[[RNA]]<-temp_Volcano
}
combine<-ggarrange(plot[["msRNA"]],plot[["lsRNA"]],plot[["miRNA"]],
                   plot[["piRNA"]],plot[["snRNA"]],plot[["snoRNA"]],
                   plot[["tsRNA"]],plot[["rsRNA"]],plot[["ysRNA"]],
                   ncol = 3, nrow = 3)
annotate_figure(combine,
                left = text_grob(expression('-Log'[10]*' adjusted p value'),rot = 90,size=18),
                bottom=text_grob(expression('Log'[2]*' fold change'),size=18))
ggsave("Volcano plot - combined 9.pdf", width=9, height=12)
ggsave("Volcano plot - combined 9.png", width=9, height=12, units="in", bg="white")

##cfRNA
results_cfRNA<-data.frame(results)
increase<-nrow(subset(results_cfRNA, padj<0.1 & log2FoldChange>0))
decrease<-nrow(subset(results_cfRNA, padj<0.1 & log2FoldChange<0))
cfRNA_Volcano<-EnhancedVolcano(results_cfRNA, lab=NA,
                               x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' adjusted p value'), 
                               title = "HCC vs. Non-HCC", subtitle = paste0("loss ",decrease,"  |  gain ",increase),caption = "",
                               pCutoff = 0.1, FCcutoff = 1, xlim = c(-9.8,6.1), ylim = c(0,26), pointSize = 1, 
                               col = c("black","black","#4E62AB","#D6404E"), colAlpha = 0.9, legendLabels = c("NS", "FC", "FDR", "FDR & FC"))
cfRNA_Volcano<-cfRNA_Volcano + theme_classic()+ theme(axis.text = element_text(color = "black"), legend.position = "none",
                                                      plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5),
                                                      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
cfRNA_Volcano
ggsave(paste0("Volcano plot - cfRNA",".pdf"), width=4, height=6)
ggsave(paste0("Volcano plot - cfRNA",".png"), width=4, height=6, units="in", bg="white")

##Distribution----
sig_results<-results_cfRNA[results_cfRNA$padj < 0.1,]
Type_number<-data.frame(group=c("msRNA","lsRNA","miRNA","piRNA","tsRNA","rsRNA","ysRNA","snRNA","snoRNA"),
                        value=c(table(rownames(sig_results) %in% rownames(msRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(lsRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(miRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(piRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(tsRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(rsRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(ysRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(snRNA))["TRUE"],
                                table(rownames(sig_results) %in% rownames(snoRNA))["TRUE"]))
Type_number[is.na(Type_number)]=0
Type_number <- Type_number %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(Type_number$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
Type_number$prop2<-paste0(round(Type_number$prop,1),"%")

Type_number$group<-factor(Type_number$group,levels=c("lsRNA","miRNA","piRNA","msRNA","snRNA","rsRNA","snoRNA","ysRNA","tsRNA"))
Type_number_out<-Type_number[Type_number$group== "miRNA" | Type_number$group== "snRNA" | Type_number$group== "snoRNA",]
Type_number_in<-Type_number[Type_number$group== "lsRNA" | Type_number$group== "piRNA" | Type_number$group== "msRNA" | Type_number$group== "rsRNA" | Type_number$group== "tsRNA" | Type_number$group== "ysRNA",]
ggplot(Type_number, aes(x="", y=prop, fill=group)) + geom_bar(stat="identity", width=1, color="black") + coord_polar("y", start=0) + 
  geom_text(data = Type_number_in, aes(y = ypos, label = paste(prop2)), color = "black", size=5, nudge_x = 0.3) +
  geom_label_repel(data = Type_number_out, aes(y = ypos, label = paste0(prop2)), size = 5, nudge_x = 0.6, show.legend = FALSE) +
  scale_fill_manual(values=c("#9E0142","#D6404E","#F57547","#FDB96A","#F0D43A","#CBE99D","#87CFA4","#469EB4","#4E62AB"))+
  theme_void() + 
  theme(plot.title=element_text(size=20,hjust = 0.5), plot.subtitle=element_text(size=16,hjust = 0.5),
        legend.text=element_text(size=14), legend.title=element_blank(), legend.position = "bottom")
ggsave(paste0("Pie plot-All cfRNA-Significant results.pdf"), width=6, height = 8)
#ggsave(paste0("Pie plot-All cfRNA-Significant results.png"), width=6, height = 8, units = "in", bg="white")
#Note: Position not correct. Edit with Adobe AI.

"#9E0142"-"rgb(158,1,66)"-"lsRNA"
"#D6404E"-"rgb(214,64,78)"-"miRNA"
"#F57547"-"rgb(245,117,71)"-"piRNA"
"#FDB96A"-"rgb(253,185,106)"-"msRNA"
"#F0D43A"-"rgb(240,212,58)"-"snRNA"
"#CBE99D"-"rgb(203,233,157)"-"rsRNA"
"#87CFA4"-"rgb(135,207,164)"-"snoRNA"
"#469EB4"-"rgb(70,158,180)"-"ysRNA"
"#4E62AB"-"rgb(78,98,171)"-"tsRNA"
"unspecify"

AllRNA<-read.table("../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[,meta$ID]
##each RNA
dir.create(paste0(getwd(),"/each RNA rpm boxplot"), showWarnings = FALSE)
subset<-AllRNA[rownames(candi_results),]
subset<-data.frame(t(subset),check.names = FALSE)
table(rownames(subset)==meta$ID)
subset$Group<-meta$Group
subset$Group<-gsub("Healthy","NOR",subset$Group)
subset$Group<-gsub("Breast cancer","BRC",subset$Group)
subset$Group<-gsub("Colorectal cancer","CRC",subset$Group)
subset$Group<-gsub("Gastric cancer","GC",subset$Group)
subset$Group<-gsub("Liver cancer","HCC",subset$Group)
subset$Group<-gsub("Lung cancer","LC",subset$Group)
subset$Group<-factor(subset$Group,levels = c("NOR","BRC","CRC","GC","HCC","LC"))
N=1:(ncol(subset)-1)
registerDoParallel(cl<-makeCluster(12))
result_seed<-foreach(i=N, .combine="c",.packages=c("tidyverse","ggplot2","ggpubr","rstatix","reshape2")) %dopar% {
  df<-subset[,c(i,ncol(subset))]
  colnames(df)[1]<-"RNA"
  stat.test <- df %>%
    t_test(RNA ~ Group, paired = FALSE, var.equal = FALSE, p.adjust.method = "none", ref.group = "HCC")%>%
    add_significance("p")%>%
    add_xy_position(fun="max")
  if (stat.test$p[1]<0.05 & stat.test$p[2]<0.05 & stat.test$p[3]<0.05 & stat.test$p[4]<0.05 & stat.test$p[5]<0.05) {
    plot <- ggplot(df, aes(x=Group, y=RNA)) + geom_boxplot(aes(fill = Group),outlier.shape=NA)+ scale_fill_manual(values=c("#4E62AB","#469EB4","#87CFA4","#F0D43A","#F57547","#D6404E"))+ 
      geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.8)+
      labs(title=paste0(colnames(subset)[i]), x="", y = expression('Log'[2]*' (RPM+1)'))+ theme_minimal()+
      theme(axis.text = element_text(colour = "black",size = 12), text = element_text(size = 14),legend.position = "none")
    plot + stat_pvalue_manual(stat.test, label = "p.signif",tip.length = 0)
    ggsave(paste0("each RNA rpm boxplot/",colnames(subset)[i],".png"),device = "png",bg="white",width = 4, height = 4)
  }}
stopCluster(cl)

  