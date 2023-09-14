rm(miRNA,piRNA,rsRNA,ysRNA,tRFs)
library(R.utils)
#library(tidyverse)
####miRNA
files = list.files(path = "result_calling/miRNA", pattern = "miRNA.txt$") 
for(file in files){
  print(file)
  File=paste("result_calling/miRNA/",file, sep="")
  if (countLines(File) <1) next
  temp = read.table(file = File)
  colnames(temp)[1] = paste0(strsplit(file, "_")[[1]][1],"_",strsplit(file, "_")[[1]][2])
  if(!exists("miRNA")){
    miRNA = temp
  } else{
    miRNA = merge(miRNA, temp, by="V2", all=T)
  }
}
row.names(miRNA)<-miRNA[,1];miRNA<-miRNA[,-1]
miRNA[is.na(miRNA)]=0
write.csv(miRNA,"result_final/miRNA.csv")

####piRNA
files = list.files(path = "result_calling/piRNA", pattern = "piRNA.txt$") 
for(file in files){
  print(file)
  File=paste("result_calling/piRNA/",file, sep="")
  if (countLines(File) <1) next
  temp = read.table(file = File)
  colnames(temp)[1] = paste0(strsplit(file, "_")[[1]][1],"_",strsplit(file, "_")[[1]][2])
  if(!exists("piRNA")){
    piRNA = temp
  } else{
    piRNA = merge(piRNA, temp, by="V2", all=T)
  }
}
row.names(piRNA)<-piRNA[,1];piRNA<-piRNA[,-1]
piRNA[is.na(piRNA)]=0
write.csv(piRNA,"result_final/piRNA.csv")

####rsRNA
files = list.files(path = "result_calling/rsRNA", pattern = "rsRNA.txt$") 
for(file in files){
  print(file)
  File=paste("result_calling/rsRNA/",file, sep="")
  if (countLines(File) <1) next
  temp = read.table(file = File)
  colnames(temp)[2] = paste0(strsplit(file, "_")[[1]][1],"_",strsplit(file, "_")[[1]][2])
  if(!exists("rsRNA")){
    rsRNA = temp
  } else{
    rsRNA = merge(rsRNA, temp, by="V1", all=T)
  }
}
row.names(rsRNA)<-rsRNA[,1];rsRNA<-rsRNA[,-1]
rsRNA[is.na(rsRNA)]=0
write.csv(rsRNA,"result_final/rsRNA.csv")

####ysRNA
files = list.files(path = "result_calling/ysRNA", pattern = "ysRNA.txt$") 
for(file in files){
  print(file)
  File=paste("result_calling/ysRNA/",file, sep="")
  if (countLines(File) <1) next
  temp = read.table(file = File)
  colnames(temp)[2] = paste0(strsplit(file, "_")[[1]][1],"_",strsplit(file, "_")[[1]][2])
  if(!exists("ysRNA")){
    ysRNA = temp
  } else{
    ysRNA = merge(ysRNA, temp, by="V1", all=T)
  }
}
row.names(ysRNA)<-ysRNA[,1];ysRNA<-ysRNA[,-1]
ysRNA[is.na(ysRNA)]=0
write.csv(ysRNA,"result_final/ysRNA.csv")

####tRFs
files = list.files(path = "result_calling/tRFs", pattern = "tRFs.txt$") 
for(file in files){
  print(file)
  File=paste("result_calling/tRFs/",file, sep="")
  if (countLines(File) <1) next
  temp = read.table(file = File)
  colnames(temp)[2] = paste0(strsplit(file, "_")[[1]][1],"_",strsplit(file, "_")[[1]][2])
  if(!exists("tRFs")){
    tRFs = temp
  } else{
    tRFs = merge(tRFs, temp, by="V1", all=T)
  }
}
row.names(tRFs)<-tRFs[,1];tRFs<-tRFs[,-1]
tRFs[is.na(tRFs)]=0
write.csv(tRFs,"result_final/tRFs.csv")

####mRNA
mRNA<-read.table("result_calling/mRNA/mRNA_S1_transcript.txt", header = T)
mRNA<-mRNA[,-(2:6)]
row.names(mRNA)<-mRNA[,1];mRNA<-mRNA[,-1]
name<-matrix(unlist(strsplit(colnames(mRNA),"\\.")),ncol=ncol(mRNA),nrow=5)
subname<-matrix(unlist(strsplit(name[4,],"_")),ncol=ncol(mRNA),nrow=3)
colnames(mRNA)<-paste0(name[3,],'-',subname[1,],"_",subname[2,])
write.csv(mRNA,"result_final/mRNA.csv")

####lncRNA
lncRNA<-read.table("result_calling/lncRNA/lncRNA_S1_transcript.txt", header = T)
lncRNA<-lncRNA[,-(2:6)]
row.names(lncRNA)<-lncRNA[,1];lncRNA<-lncRNA[,-1]
name<-matrix(unlist(strsplit(colnames(lncRNA),"\\.")),ncol=ncol(lncRNA),nrow=5)
subname<-matrix(unlist(strsplit(name[4,],"_")),ncol=ncol(lncRNA),nrow=3)
colnames(lncRNA)<-paste0(name[3,],'-',subname[1,],"_",subname[2,])
write.csv(lncRNA,"result_final/lncRNA.csv")

####snRNA
snRNA<-read.table("result_calling/snRNA/snRNA_S1_transcript.txt", header = T)
snRNA<-snRNA[,-(2:6)]
row.names(snRNA)<-snRNA[,1];snRNA<-snRNA[,-1]
name<-matrix(unlist(strsplit(colnames(snRNA),"\\.")),ncol=ncol(snRNA),nrow=5)
subname<-matrix(unlist(strsplit(name[4,],"_")),ncol=ncol(snRNA),nrow=3)
colnames(snRNA)<-paste0(name[3,],'-',subname[1,],"_",subname[2,])
write.csv(snRNA,"result_final/snRNA.csv")

####snoRNA
snoRNA<-read.table("result_calling/snoRNA/snoRNA_S1_transcript.txt", header = T)
snoRNA<-snoRNA[,-(2:6)]
row.names(snoRNA)<-snoRNA[,1];snoRNA<-snoRNA[,-1]
name<-matrix(unlist(strsplit(colnames(snoRNA),"\\.")),ncol=ncol(snoRNA),nrow=5)
subname<-matrix(unlist(strsplit(name[4,],"_")),ncol=ncol(snoRNA),nrow=3)
colnames(snoRNA)<-paste0(name[3,],'-',subname[1,],"_",subname[2,])
write.csv(snoRNA,"result_final/snoRNA.csv")
