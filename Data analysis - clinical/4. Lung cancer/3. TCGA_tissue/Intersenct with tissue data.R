library(tidyverse)
library(ggvenn)
library(ggpubr)
library(miRBaseConverter)

RNA_name <- list.files(path = "../../Normalization_RPM_colsum/", pattern = "cpm.txt$", full.names = TRUE, recursive = TRUE) 
RNA_name <- RNA_name[grep("cfRNA",RNA_name, invert = TRUE)]
RNA_name <- lapply(RNA_name, read.table, header=TRUE, sep="\t",row.names=1)
RNA_name <- lapply(RNA_name, row.names)
names(RNA_name)<- c("lsRNA",  "miRNA",  "msRNA",  "piRNA",  "rsRNA",  "snoRNA", "snRNA",  "tsRNA",  "ysRNA")
sig_cfRNA<-read.table("../DE_LC_SZDE_vs_NOR_SZBA/DESeq2 significant results-cfRNA-Lung cancer vs. Healthy.txt")
sig_list<-list()
Distribution<-double()
for (RNA in names(RNA_name)) {
  sig_list[[RNA]]<-rownames(sig_cfRNA)[(rownames(sig_cfRNA) %in% RNA_name[[RNA]])]
  Distribution[RNA]=table(rownames(sig_cfRNA) %in% RNA_name[[RNA]])["TRUE"]
}
sum(Distribution)


####mRNA----
plasma_sig<-sig_list[["msRNA"]]
tissue_sig<-rownames(read.table("DESeq2 significant results-TCGA-LUAD-mRNA-Primary Tumor vs. Solid Tissue Normal.txt"))
table(plasma_sig %in% tissue_sig)
# FALSE  TRUE 
# 454  1023 
table(tissue_sig %in% plasma_sig)
# FALSE  TRUE 
# 12361  1023 
Venn_list<-list(
  Plasma=plasma_sig,
  Tissue=tissue_sig)
names(Venn_list) <- c("Plasma msRNA","Tissue mRNA")
Venn <- ggvenn(Venn_list, show_percentage = F, set_name_size = 8, text_size = 7,fill_color = c("red", "yellow"),fill_alpha = 0.5)
ggsave("Venn diagram-Plasma msRNA & Tissue mRNA.pdf", width=6, height = 6)
ggsave("Venn diagram-Plasma msRNA & Tissue mRNA.png", width=6, height = 6, bg = "white")


####miRNA----
plasma_sig<-sig_list[["miRNA"]]
plasma_version=checkMiRNAVersion(plasma_sig, verbose = TRUE)
plasma_new<-miRNA_MatureToPrecursor(plasma_sig, version=plasma_version)
tissue_sig<-rownames(read.table("DESeq2 significant results-TCGA-LUAD-miRNA-Primary Tumor vs. Solid Tissue Normal.txt"))
tissue_version=checkMiRNAVersion(tissue_sig, verbose = TRUE)
tissue_new<-miRNAVersionConvert(tissue_sig,targetVersion = plasma_version,exact = TRUE)
table(plasma_new$Precursor %in% tissue_new$TargetName)
# FALSE  TRUE 
# 140   173  
#table(tissue_new$TargetName %in% plasma_new$Precursor)
# FALSE  TRUE 
# 335   134
Venn_list<-list(
  Plasma = plasma_new$Precursor,
  Tissue = tissue_new$TargetName)
names(Venn_list) <- c("Plasma miRNA","Tissue miRNA")
Venn <- ggvenn(Venn_list, show_percentage = F, set_name_size = 8, text_size = 7,fill_color = c("red", "yellow"),fill_alpha = 0.5)
ggsave("Venn diagram-Plasma miRNA & Tissue miRNA.pdf", width=6, height = 6)
#Edit with Adobe Illustrator for incorrect numbers due to duplicate precursor names
