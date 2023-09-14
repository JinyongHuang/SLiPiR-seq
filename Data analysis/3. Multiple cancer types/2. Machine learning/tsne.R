library(tidyverse)
library(ggplot2)
library(ggVennDiagram)
library(Rtsne)
library(rgl)
####Sample ID----
meta<-read.csv("../../Metadata/Meta_with_clinical_and_Summary_good_samples_20230829.csv", header = T)
meta <- subset(meta,sample_type=="LC_SZDE" | 
                 sample_type=="BRC_NXYK" | sample_type=="CRC_NXYK" | sample_type=="GC_NXYK" | sample_type=="HCC_NXYK")

###GOI----
BRC<-readRDS("../ML_BRC_vs._others/goi_list_Breast cancer vs. Non-Breast cancer.rds")
CRC<-readRDS("../ML_CRC_vs._others/goi_list_Colorectal cancer vs. Non-Colorectal cancer.rds")
GC<-readRDS("../ML_GC_vs._others/goi_list_Gastric cancer vs. Non-Gastric cancer.rds")
HCC<-readRDS("../ML_HCC_vs._others/goi_list_Liver cancer vs. Non-Liver cancer.rds")
LC<-readRDS("../ML_LC_vs._others/goi_list_Lung cancer vs. Non-Lung cancer.rds")
goi_list<-list(BRC=BRC,CRC=CRC,GC=GC,HCC=HCC,LC=LC)
goi_list<-lapply(goi_list, "[[", "cfRNA")
goi<-c(goi_list[[1]],goi_list[[2]],goi_list[[3]],goi_list[[4]],goi_list[[5]])
goi<-unique(goi)

###data----
AllRNA<-read.table("../../Normalization_RPM_colsum/cfRNA-log2(cpm).txt", header = T, row.names = 1)
AllRNA<-AllRNA[goi,meta$ID]

####2D 2D tsne
tsne_df<-data.frame(t(AllRNA),check.names = FALSE)
for (i in 1:10){
  set.seed(i)
  tsne <- Rtsne(tsne_df, dims = 2, perplexity=10, check_duplicates = FALSE)
  tsne <-as.data.frame(tsne$Y)
  rownames(tsne)<-meta$ID
  colnames(tsne)<-c("tSNE1","tSNE2")
  table(rownames(tsne)==meta$ID)
  tsne$Group<-factor(meta$Group)
ggplot(tsne, aes(x=tSNE1, y=tSNE2,col= Group,fill=Group)) + 
  geom_point() + theme_classic()+ labs(x="t-SNE 1",y="t-SNE 2")+
  stat_ellipse(level=0.6, geom="polygon", alpha=1/3)+
  scale_fill_manual(values =c("#469EB4","#87CFA4","#F0D43A","#F57547","#D6404E"))+
  scale_color_manual(values =c("#469EB4","#87CFA4","#F0D43A","#F57547","#D6404E"))+
  theme(legend.title = element_blank(),legend.position = c(0.15,0.15),legend.background = element_blank(),
        axis.text = element_text(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3),
        axis.line = element_line(linewidth=0.15), text = element_text(size=10))
ggsave(paste0("tsne/2D tsne_",i,".png"),width = 5, height = 5, bg="white")
ggsave(paste0("tsne/2D tsne_",i,".pdf"),width = 5, height = 5)
}
####3D tsne
# meta$Color<-NA
# meta$Color[which(meta$Group=="Breast cancer")]<-"#469EB4"
# meta$Color[which(meta$Group=="Colorectal cancer")]<-"#87CFA4"
# meta$Color[which(meta$Group=="Gastric cancer")]<-"#F0D43A"
# meta$Color[which(meta$Group=="Liver cancer")]<-"#F57547"
# meta$Color[which(meta$Group=="Lung cancer")]<-"#D6404E"
# colors <- meta$Color
# tsne <- Rtsne(tsne_df, dims = 3, perplexity=20, check_duplicates = FALSE)
# plot<-plot3d(x=tsne$Y[,1], y=tsne$Y[,2], z=tsne$Y[,3], type = 'p', col = colors,size=6)
