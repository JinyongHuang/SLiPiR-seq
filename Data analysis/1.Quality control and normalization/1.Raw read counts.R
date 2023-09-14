library(tidyverse)
library(DescTools)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(rstatix)

####Total----
rm(Total)
files = list.files(path = "../Raw data/", pattern = "Total_reads.txt", full.names = TRUE, recursive = TRUE) 
for (file in files) {
  print(file)
  temp<-read.table(file)
  name<-matrix(unlist(strsplit(temp$V2,"/")),ncol=nrow(temp),nrow=3)[3,]
  name<-matrix(unlist(strsplit(name,"_")),ncol=nrow(temp),nrow=4)
  row.names(temp)<-paste0(name[1,],"_",name[2,])
  temp<-temp[1]
  colnames(temp)<-"Total"
  if (exists("Total")){Total<-rbind(Total,temp)}
  if (!exists("Total")){Total<-temp}
  rm(temp)
}
Total<-Total[OrderMixed(row.names(Total)), , drop=FALSE]

####Trimmed
rm(Trimmed)
files = list.files(path = "../Raw data/", pattern = "Trimmed_reads.txt", full.names = TRUE, recursive = TRUE) 
for (file in files) {
  print(file)
  temp<-read.table(file)
  name<-matrix(unlist(strsplit(temp$V2,"/")),ncol=nrow(temp),nrow=3)[3,]
  name<-matrix(unlist(strsplit(name,"_")),ncol=nrow(temp),nrow=5)
  row.names(temp)<-paste0(name[1,],"_",name[2,])
  temp<-temp[order(row.names(temp)), ]
  temp<-temp[1]
  colnames(temp)<-"Trimmed"
  if (exists("Trimmed")){Trimmed<-rbind(Trimmed,temp)}
  if (!exists("Trimmed")){Trimmed<-temp}
  rm(temp)
}
Trimmed<-Trimmed[OrderMixed(row.names(Trimmed)), , drop=FALSE]

####Unmapped
rm(Unmapped)
files = list.files(path = "../Raw data/", pattern = "Unmapped_reads.txt", full.names = TRUE, recursive = TRUE) 
for (file in files) {
  print(file)
  temp<-read.table(file)
  name<-matrix(unlist(strsplit(temp$V2,"/")),ncol=nrow(temp),nrow=3)[3,]
  name<-matrix(unlist(strsplit(name,"_")),ncol=nrow(temp),nrow=3)
  row.names(temp)<-paste0(name[1,],"_",name[2,])
  temp<-temp[order(row.names(temp)), ]
  temp<-temp[1]
  colnames(temp)<-"Unmapped"
  if (exists("Unmapped")){Unmapped<-rbind(Unmapped,temp)}
  if (!exists("Unmapped")){Unmapped<-temp}
  rm(temp)
}
Unmapped<-Unmapped[OrderMixed(row.names(Unmapped)), , drop=FALSE]

Summary<-cbind(Total,Trimmed,Unmapped)
Summary$CleanRatio<-round(Summary$Trimmed/Summary$Total*100,2)
Summary$MapRate<-round((Summary$Total-Summary$Unmapped)/Summary$Total*100,2)

####RNA read counts----
RNAs=c("lncRNA", "miRNA", "mRNA", "piRNA", "rsRNA", "snoRNA", "snRNA", "tRFs", "ysRNA" )
rawdata<-list()
for (RNA in RNAs) {
  files = list.files(path = "../Raw data/", pattern = paste0(RNA,".csv"), full.names = TRUE, recursive = TRUE)
  for (file in files) {
    print(file)
    if (exists("temp1")){
      temp2 <- read.csv(file = file, row.names = 1, check.names = FALSE)
      temp1 <- merge(temp1,temp2,by=0,all=T)
      rownames(temp1)<-temp1[,1]
      temp1<-temp1[,-1]
    }
    if (!exists("temp1")){
      temp1<-read.csv(file = file, row.names = 1, check.names = FALSE)
    }
  }
  temp1[is.na(temp1)] <- 0 
  rawdata[[RNA]]=temp1
  rm(temp1,temp2)
}


####RNA sum and ratio
#Sum<-lapply(rawdata,function(x) apply(x, 2, sum))
#Sum<-as.data.frame(do.call(cbind, Sum))
rm(Sum)
for (RNA in RNAs) {
  temp<-rawdata[[RNA]]
  temp_sum<-data.frame(apply(temp, 2, sum))
  colnames(temp_sum)<-paste0(RNA,"_sum")
  if (exists("Sum")){
    Sum=merge(Sum,temp_sum,by=0,all=T)
    rownames(Sum)<-Sum[,1]
    Sum<-Sum[,-1]
  }
  if (!exists("Sum")){
    Sum=temp_sum
  }
}
Sum[is.na(Sum)] <- 0 
Sum<-Sum[OrderMixed(rownames(Sum)),]
Sum$cfRNA_sum<-apply(Sum,1,sum)
Summary<-cbind(Summary,Sum)
ratio<-round(Summary[,grep("sum",colnames(Summary))]/Summary$Trimmed*100,2)
colnames(ratio)<-gsub("sum","ratio",colnames(ratio))
ratio$Unspecified_ratio<-100-ratio$cfRNA_ratio
Summary<-cbind(Summary,ratio)

####RNA detected
rm(detected)
for (RNA in RNAs) {
  temp<-rawdata[[RNA]]
  temp_detected<-data.frame(apply(temp, 2, function(x) length(which(x!=0))))
  colnames(temp_detected)<-paste0(RNA,"_detected")
  if (exists("detected")){
    detected=merge(detected,temp_detected,by=0,all=T)
    rownames(detected)<-detected[,1]
    detected<-detected[,-1]
  }
  if (!exists("detected")){
    detected=temp_detected
  } 
}
detected[is.na(detected)] <- 0 
detected<-detected[OrderMixed(rownames(detected)),]
detected$cfRNA_detected<-apply(detected,1,sum)
Summary<-cbind(Summary,detected)

write.csv(Summary,"Summary_sequencing_submission.csv")

####Unique samples----
meta_ID<-read.csv("../Metadata/Meta_id_seqID.csv", header = T)
meta_ID<-meta_ID[OrderMixed(meta_ID$id),]
Summary$seq_ID<-meta_ID$seq_ID
uni_summary<-aggregate(Summary[,1:3],by=list(Summary$seq_ID),sum)
rownames(uni_summary)<-uni_summary[,1];uni_summary<-uni_summary[,-1]
uni_summary$CleanRatio<-round(uni_summary$Trimmed/uni_summary$Total*100,2)
uni_summary$MapRate<-round((uni_summary$Total-uni_summary$Unmapped)/uni_summary$Total*100,2)
uni_summary<-uni_summary[OrderMixed(rownames(uni_summary)),]
uni_rawdata<-list()
for (RNA in RNAs) {
  print(RNA)
  temp<-rawdata[[RNA]]
  temp<-data.frame(t(temp),check.names = FALSE)
  temp<-temp[OrderMixed(rownames(temp)),]
  temp$seq_ID<-meta_ID$seq_ID
  temp1<-aggregate(temp[,-ncol(temp)],by=list(temp$seq_ID),sum)
  rownames(temp1)<-temp1[,1];temp1<-temp1[,-1]
  temp1<-data.frame(t(temp1),check.names = FALSE)
  uni_rawdata[[RNA]]<-temp1
  rm(temp,temp1)
}
####RNA sum and ratio
rm(uni_Sum)
for (RNA in RNAs) {
  temp<-uni_rawdata[[RNA]]
  temp_sum<-data.frame(apply(temp, 2, sum))
  colnames(temp_sum)<-paste0(RNA,"_sum")
  if (exists("uni_Sum")){
    uni_Sum=merge(uni_Sum,temp_sum,by=0,all=T)
    rownames(uni_Sum)<-uni_Sum[,1]
    uni_Sum<-uni_Sum[,-1]
  }
  if (!exists("uni_Sum")){
    uni_Sum=temp_sum
  }
}
uni_Sum[is.na(uni_Sum)] <- 0 
uni_Sum<-uni_Sum[OrderMixed(rownames(uni_Sum)),]
uni_Sum$cfRNA_sum<-apply(uni_Sum,1,sum)
uni_summary<-cbind(uni_summary,uni_Sum)
uni_ratio<-round(uni_summary[,grep("sum",colnames(uni_summary))]/uni_summary$Trimmed*100,2)
colnames(uni_ratio)<-gsub("sum","ratio",colnames(uni_ratio))
uni_ratio$Unspecified_ratio<-100-uni_ratio$cfRNA_ratio
uni_summary<-cbind(uni_summary,uni_ratio)

####RNA detected
rm(uni_detected)
for (RNA in RNAs) {
  temp<-uni_rawdata[[RNA]]
  temp_detected<-data.frame(apply(temp, 2, function(x) length(which(x!=0))))
  colnames(temp_detected)<-paste0(RNA,"_detected")
  if (exists("uni_detected")){
    uni_detected=merge(uni_detected,temp_detected,by=0,all=T)
    rownames(uni_detected)<-uni_detected[,1]
    uni_detected<-uni_detected[,-1]
  }
  if (!exists("uni_detected")){
    uni_detected=temp_detected
  } 
}
uni_detected[is.na(uni_detected)] <- 0 
uni_detected<-uni_detected[OrderMixed(rownames(uni_detected)),]
uni_detected$cfRNA_detected<-apply(uni_detected,1,sum)
uni_summary<-cbind(uni_summary,uni_detected)

write.csv(uni_summary,"Summary_unique_sample.csv")

meta_unique<-read.csv("../Metadata/Meta_clinical.csv", header = T)
meta_unique<-meta_unique[OrderMixed(meta_unique$sample_type),]
ntype<-data.frame(table(meta_unique$sample_type))
name<-apply(ntype, 1, function(x) str_c(x[1],"_",1:x[2]))
for (i in 1:length(name)) {
  if (i==1) {ID=name[[i]]}
  if (i!=1) {ID=c(ID,name[[i]])}
}
meta_unique$ID<-ID

meta<-merge(meta_unique,uni_summary,by.x="seq_ID",by.y = 0)
meta<-meta[OrderMixed(meta$seq_ID),]
write.csv(meta,"../Metadata/Meta_with_clinical_and_Summary_all_samples.csv",row.names = F)

####convert sample ID----
for (RNA in RNAs) {
  print(RNA)
  temp<-uni_rawdata[[RNA]]
  temp<-temp[,OrderMixed(colnames(temp))]
  table(meta$seq_ID==colnames(temp))
  colnames(temp)=meta$ID
  temp<-temp[,OrderMixed(colnames(temp))]
  write.csv(temp,paste0(RNA,".csv"))
}

####Filter out bad samples
meta_good<-meta_unique[meta_unique$CleanRatio>20,]#bad library
meta_good<-meta_good[meta_good$Trimmed>2e6,]#low informative reads
meta_good<-meta_good[meta_good$rsRNA_ratio<30,]#rsRNA contemination
meta_good$gRNA_ratio<-meta_good$lncRNA_ratio+meta_good$mRNA_ratio
meta_good<-meta_good[meta_good$gRNA_ratio<30,]#cell contemination
meta_good<-meta_good[-(which(meta_good$sample_type=="LC_SZBU" & meta_good$AJCC.Stage=="IV")),]#stage IV in validation cohort
table(meta_good$sample_type)
write.csv(meta_good,"../Metadata/Meta_with_clinical_and_Summary_good_samples.csv",row.names = F)

####Used samples----
meta_use<-meta_good[meta_good$sample_type=="BRC_NXYK" | meta_good$sample_type=="CRC_NXYK" | meta_good$sample_type=="GC_NXYK" | meta_good$sample_type=="HCC_NXYK" |
                      meta_good$sample_type=="LC_SZDE" | meta_good$sample_type=="NOR_SZBA" | meta_good$sample_type=="LC_SZBU" | meta_good$sample_type=="NOR_SZDW",]
write.csv(meta_use,"../Metadata/Meta_with_clinical_and_Summary_used_samples.csv",row.names = F)

##raw reads of used samples for publication
dir.create(paste0(getwd(),"/Data submission"), showWarnings = FALSE)
meta_good<-read.csv("../Metadata/Meta_with_clinical_and_Summary_good_samples.csv")
meta_use<-meta_good[meta_good$sample_type=="BRC_NXYK" | meta_good$sample_type=="CRC_NXYK" | meta_good$sample_type=="GC_NXYK" | meta_good$sample_type=="HCC_NXYK" |
                      meta_good$sample_type=="LC_SZDE" | meta_good$sample_type=="NOR_SZBA" | meta_good$sample_type=="LC_SZBU" | meta_good$sample_type=="NOR_SZDW",]
meta_use<-meta_use[OrderMixed(meta_use$seq_ID),]

for (RNA in RNAs) {
  temp<-uni_rawdata[[RNA]]
  temp<-temp[,meta_use$seq_ID]
  temp<-temp[,OrderMixed(colnames(temp))]
  table(meta_use$seq_ID==colnames(temp))
  colnames(temp)=meta_use$ID
  temp<-temp[,OrderMixed(colnames(temp))]
  write.csv(temp,paste0("Data submission/",RNA,".csv"))
}
