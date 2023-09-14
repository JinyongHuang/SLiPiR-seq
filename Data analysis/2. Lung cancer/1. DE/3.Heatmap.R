library(tidyverse)
library(DescTools)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(circlize)#colorRamp2
library(RColorBrewer)
library(ComplexHeatmap)#Heatmap

####Sample ID----
meta<-read.csv("../../Metadata/Meta_with_clinical_and_Summary_good_samples_20230829.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="NOR_SZBA")

####Significant RNA
RNA_name <- list.files(path = "../../Normalization_RPM_colsum/", pattern = "cpm.txt$", full.names = TRUE, recursive = TRUE) 
RNA_name <- RNA_name[grep("cfRNA",RNA_name, invert = TRUE)]
RNA_name <- lapply(RNA_name, read.table, header=TRUE, sep="\t",row.names=1)
RNA_name <- lapply(RNA_name, row.names)
names(RNA_name)<- c("lncRNA",  "miRNA",  "mRNA",  "piRNA",  "rsRNA",  "snoRNA", "snRNA",  "tsRNA",  "ysRNA")
sig_cfRNA<-read.table("../DE_LC_SZDE_vs_NOR_SZBA/DESeq2 candidate results-cfRNA-Lung cancer vs. Healthy.txt", header = T)
sig_cfRNA<-sig_cfRNA[sig_cfRNA$log2FoldChange>=0,]
sig_list<-list()
sig_name_list<-list()
top_name_list<-list()
Distribution<-double()
for (RNA in names(RNA_name)) {
  sig_name<-rownames(sig_cfRNA)[(rownames(sig_cfRNA) %in% RNA_name[[RNA]])]
  sig_name_list[[RNA]]<-sig_name
  sig_list[[RNA]]<-sig_cfRNA[sig_name,]
  top_name_list[[RNA]]<-rownames(sig_cfRNA[sig_name,][1:15,])
  Distribution[RNA]=table(rownames(sig_cfRNA) %in% RNA_name[[RNA]])["TRUE"]
}
GOI<-unlist(top_name_list, use.names = FALSE)
###RNA reads data
AllRNA<-read.table("../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T)
AllRNA<-AllRNA[GOI,meta$ID]
meta<-meta[OrderMixed(meta$ID),]
AllRNA<-AllRNA[,OrderMixed(colnames(AllRNA))]

#heatmap-log2RPM
# df<-as.matrix(AllRNA)
# type <- meta$Group
# col_type <- c("Healthy" = "green", "Lung cancer" = "orange")
# hA <- HeatmapAnnotation(Type = type, col = list(Type = col_type), show_annotation_name = F)
# h <- Heatmap(df, cluster_rows = F, show_row_dend = T, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
#              cluster_columns = F, show_column_dend = T, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
#              col = colorRamp2(c(0, max(df)), c("white", "red")), show_column_names = FALSE,
#              column_labels = colnames(df), column_title = "Hierarchical clustering (Euclidean distance, Ward's method)", 
#              name="log2RPM",top_annotation = hA)
# print(h)

##heatmap-z.score
# z.score <- t(apply(AllRNA,1,scale)) #z-score
# colnames(z.score) <- colnames(AllRNA)
# type <- meta$Group
# col_type <- c("Healthy" = "green", "Lung cancer" = "orange")
# top<- HeatmapAnnotation(Type = type, col = list(Type = col_type), show_annotation_name = F)
# h <- Heatmap(z.score, cluster_rows = T, show_row_dend = T, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
#              cluster_columns = T, show_column_dend = T, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
#              col = colorRamp2(c(quantile(z.score,0), 0, quantile(z.score,1)), c("blue", "white", "red")), show_column_names = FALSE,
#              column_title = "Hierarchical clustering (Euclidean distance, Ward's method)", 
#              name="z-score",top_annotation = top)
# print(h)

             
##heatmap-z.score
z.score <- t(apply(AllRNA,1,scale)) #z-score
colnames(z.score) <- colnames(AllRNA)
type <- meta$Group
col_type <- c("Healthy" = "#4E62AB", "Lung cancer" = "#D6404E")
top<- HeatmapAnnotation(Group = type, col = list(Group = col_type), show_annotation_name = F, which="column")
RNA_type<-c(rep("lncRNA",15),rep("miRNA",15),rep("mRNA",15),rep("piRNA",15),
            rep("rsRNA",15),rep("snoRNA",15),rep("snRNA",15),rep("tsRNA",15),rep("ysRNA",15))
RNA_col<-c("lncRNA"="#9E0142","miRNA"="#D6404E","piRNA"="#F57547","mRNA"="#FDB96A",
           "snRNA"="#F0D43A","rsRNA"="#CBE99D","snoRNA"="#87CFA4","ysRNA"="#469EB4","tsRNA"="#4E62AB")
left<-HeatmapAnnotation(RNA = RNA_type, col = list(RNA = RNA_col), show_annotation_name = F, which="row")
for (i in 1:nrow(meta)) {
  if (meta$Group[i]=="Lung cancer" & is.na(meta$AJCC.Stage[i])) {meta$AJCC.Stage[i]="missing"}
  if (meta$Group[i]=="Healthy" & is.na(meta$AJCC.Stage[i])) {meta$AJCC.Stage[i]="NOR"}
  
}
Stage_type<-meta$AJCC.Stage
Stage_col<-c("I"="#FB7B5B","II"="#F5553C","III"="#E32F27","IV"="#C2161B","missing"="#BDBDBD","NOR"="white")
bottom<-HeatmapAnnotation(Stage = Stage_type, col = list(Stage = Stage_col), show_annotation_name = F, which="column")

h<-Heatmap(z.score, cluster_rows = T, show_row_dend = F, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
              cluster_columns = F, show_column_dend = T, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
              col = colorRamp2(c(quantile(z.score,0.01), 0, quantile(z.score,0.99)), c("#4E62AB", "white", "#9E0142")), 
             show_column_names = FALSE, show_row_names = FALSE, 
             column_title = NULL, 
             name="z-score",top_annotation = top, left_annotation =left, bottom_annotation  = bottom)
print(h)
# Export with R
# png(file="Hierarchical clustering (Euclidean distance, Ward's method).png",res=300, width=ncol(z.score)*20, height=nrow(z.score)*60)
# pdf(file="Hierarchical clustering (Euclidean distance, Ward's method).pdf",width=ncol(z.score)*0.3, height=nrow(z.score)*0.1)
# print(h);dev.off()
