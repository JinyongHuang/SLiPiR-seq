library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(RColorBrewer)
library(circlize)#colorRamp2
library(patchwork)#Combine ggplot
library(DESeq2)#Differential analysis
library(pcaExplorer)#PCA
library(EnhancedVolcano)#Volcano

####metadata----
meta<-read.csv("gdc_sample_sheet.2022-09-30.tsv",header = T, sep = "\t")
table(duplicated(meta$Case.ID))
STN<-subset(meta,Sample.Type=="Solid Tissue Normal")
PT<-subset(meta,Sample.Type=="Primary Tumor")
PT<-PT[PT$Case.ID %in% STN$Case.ID,]
PT<-PT[grepl("01A", PT$Sample.ID),]
PT<-PT[!duplicated(PT$Case.ID),]
STN<-STN[STN$Case.ID %in% PT$Case.ID,]
STN<-STN[!duplicated(STN$Case.ID),]
paired<-rbind(PT,STN)
paired$Sample.Type<-factor(paired$Sample.Type,levels = c("Solid Tissue Normal","Primary Tumor"))

####read count----
files = list.files(path = "data/.", pattern = "mirbase21.mirnas.quantification.txt$", recursive = TRUE) 
for(file in files){
  print(file)
  temp = fread(file = paste("data/",file, sep=""))
  temp = temp[,1:2]
  colnames(temp)[2] = strsplit(file, "/")[[1]][2]
  if(!exists("tcga_LUAD")){
    tcga_LUAD = temp
  } else{
    tcga_LUAD = merge(tcga_LUAD, temp, by="miRNA_ID")
  }
}
rm(temp)
tcga_LUAD = as.data.frame(tcga_LUAD)
row.names(tcga_LUAD)<-tcga_LUAD[,1]
tcga_LUAD<-tcga_LUAD[,-1]

tcga_LUAD_paired<-tcga_LUAD[,paired$File.Name]
colnames(tcga_LUAD_paired)[which(colnames(tcga_LUAD_paired) %in% PT$File.Name)]<-str_c("PT_",1:nrow(PT))
colnames(tcga_LUAD_paired)[which(colnames(tcga_LUAD_paired) %in% STN$File.Name)]<-str_c("STN_",1:nrow(STN))
write.table(tcga_LUAD_paired,"TCGA_LUAD_miRNA_raw_reads.txt", quote = F,  sep = "\t", col.names = NA)

####DE----
case="Primary Tumor"
control="Solid Tissue Normal"
ncase=table(paired$Sample.Type==case)["TRUE"]
ncontrol=table(paired$Sample.Type==control)["TRUE"]

dds<-DESeqDataSetFromMatrix(tcga_LUAD_paired, paired, ~ Sample.Type)
dds<-DESeq(dds)
results<-results(dds, contrast = c("Sample.Type", case, control)) 
summary(results)
results_df<-data.frame(results)
MA<-ggmaplot(results_df, fdr=0.1, fc=(2^0.5), size=1, top = 0, legend = "top",
             main="TCGA-LUAD-miRNA-DESeq2 results", 
             xlab = expression('Log'[2]*' mean expression'), ylab=expression('Log'[2]*' fold change'), ggtheme = ggplot2::theme_minimal())
MA<-MA+labs(subtitle = paste0(case," (N=",ncase,") vs. ",control," (N=",ncontrol,")"))
png(file=paste0("MA plot-TCGA-LUAD-miRNA-",case," vs. ",control,".png"),res=150, width=800, height=800)
print(MA);dev.off()
increase<-nrow(subset(data.frame(results_df), padj<0.1 & log2FoldChange>0))
decrease<-nrow(subset(data.frame(results_df), padj<0.1 & log2FoldChange<0))
Volcano<-EnhancedVolcano(results_df, lab = NA, x = "log2FoldChange", y = "padj", ylab = expression('-Log'[10]*' adjusted p value'), 
                         title = "TCGA-LUAD-miRNA", subtitle = paste0("loss ",decrease,"  |  gain ",increase),caption=NULL,
                         pCutoff = 0.1, FCcutoff = 1, xlim = c(-6.6,7.4), ylim = c(0,48), 
                         pointSize = 0.3, colAlpha = 0.9, col = c("black","black","#4E62AB","#D6404E"))
Volcano<-Volcano + theme_classic() + theme(axis.text = element_text(color = "black",size=6), axis.title = element_text(size=6),
                        plot.title = element_text(size=6,hjust=0.5,face="bold"), plot.subtitle = element_text(size=6,hjust=0.5),
                        legend.position = "none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3))
Volcano
ggsave(paste0("Volcano plot-TCGA-LUAD-miRNA-",case," vs. ",control,".pdf"), width=4, height=6, units="cm")
ggsave(paste0("Volcano plot-TCGA-LUAD-miRNA-",case," vs. ",control,".png"), width=4, height=6, units="cm", bg="white")


results_df<-results_df[!(is.na(results_df$padj)),]
results_df<-results_df[order(results_df$pvalue),]
write.table(results_df, paste0("DESeq2 results-TCGA-LUAD-miRNA-",case," vs. ",control,".txt"),sep = "\t", col.names = NA, quote = FALSE)
sig_results<-results_df[results_df$padj < 0.1,]
write.table(sig_results, paste0("DESeq2 significant results-TCGA-LUAD-miRNA-",case," vs. ",control,".txt"),sep = "\t", col.names = NA, quote = FALSE)
