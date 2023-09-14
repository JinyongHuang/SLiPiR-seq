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
meta <- subset(meta,sample_type=="LC_SZDE" | sample_type=="BRC_NXYK" | sample_type=="CRC_NXYK" | sample_type=="GC_NXYK" | sample_type=="HCC_NXYK")

###GOI----
BRC<-readRDS("../ML_BRC_vs._others/goi_list_Breast cancer vs. Non-Breast cancer.rds")
CRC<-readRDS("../ML_CRC_vs._others/goi_list_Colorectal cancer vs. Non-Colorectal cancer.rds")
GC<-readRDS("../ML_GC_vs._others/goi_list_Gastric cancer vs. Non-Gastric cancer.rds")
HCC<-readRDS("../ML_HCC_vs._others/goi_list_Liver cancer vs. Non-Liver cancer.rds")
LC<-readRDS("../ML_LC_vs._others/goi_list_Lung cancer vs. Non-Lung cancer.rds")
goi_list<-list(BRC=BRC,CRC=CRC,GC=GC,HCC=HCC,LC=LC)
goi_list<-lapply(goi_list, "[[", "cfRNA")
goi<-c(goi_list[[1]],goi_list[[2]],goi_list[[3]],goi_list[[4]],goi_list[[5]])
goi<-names(which(table(goi)==1))

nBRC<-table(goi %in% goi_list[[1]])["TRUE"]
goiBRC<-goi[goi %in% goi_list[[1]]]
nCRC<-table(goi %in% goi_list[[2]])["TRUE"]
goiCRC<-goi[goi %in% goi_list[[2]]]
nGC<-table(goi %in% goi_list[[3]])["TRUE"]
goiGC<-goi[goi %in% goi_list[[3]]]
nHCC<-table(goi %in% goi_list[[4]])["TRUE"]
goiHCC<-goi[goi %in% goi_list[[4]]]
nLC<-table(goi %in% goi_list[[5]])["TRUE"]
goiLC<-goi[goi %in% goi_list[[5]]]
nBRC+nCRC+nGC+nHCC+nLC
GOI<-c(goiBRC,goiCRC,goiGC,goiHCC,goiLC)
###RNA reads data
AllRNA<-read.table("../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T)
AllRNA<-AllRNA[GOI,meta$ID]
meta<-meta[OrderMixed(meta$ID),]
AllRNA<-AllRNA[,OrderMixed(colnames(AllRNA))]

           
##heatmap-z.score
z.score <- t(apply(AllRNA,1,scale)) #z-score
colnames(z.score) <- colnames(AllRNA)
type <- meta$Group
col_type <- c("Breast cancer"="#469EB4", "Colorectal cancer"="#87CFA4", "Gastric cancer"="#F0D43A", "Liver cancer"="#F57547", "Lung cancer"="#D6404E")
top<- HeatmapAnnotation(Group = type, col = list(Group = col_type), show_annotation_name = F, which="column")
panel_type<-c(rep("BRC panel",nBRC),rep("CRC panel",nCRC),rep("GC panel",nGC),rep("HCC panel",nHCC),rep("LC panel",nLC))
panel_col <- c("BRC panel"="#469EB4", "CRC panel"="#87CFA4", "GC panel"="#F0D43A", "HCC panel"="#F57547", "LC panel"="#D6404E")
left<-HeatmapAnnotation(Panel = panel_type, col = list(Panel = panel_col), show_annotation_name = F, which="row")

#methods<-c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
h<-Heatmap(z.score, cluster_rows = F, cluster_columns = F, 
              col = colorRamp2(c(quantile(z.score,0.01), 0, quantile(z.score,0.99)), c("#4E62AB", "white", "#9E0142")), 
             show_column_names = FALSE, show_row_names = FALSE, 
             column_title = NULL, name="z-score",top_annotation = top, left_annotation =left)
print(h)
#bad
# Export with R
# show_row_dend = F, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
# show_column_dend = T, clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
