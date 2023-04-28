library(tidyverse)
library(microseq)
library(foreach)#multi-core to save time
library(doParallel)#multi-core to save time

####tRFs----
tRFs<-read.table("../Reference/tRFs/tRFs.txt", sep = "\t")#https://cm.jefferson.edu/MINTbase/
colnames(tRFs)<-c("License Plate","Fragment sequence","Fragment Length","5'-half","5'-tRF","i-tRF","3'-tRF","3'-half",
                       "Total number of genomic instances in known tRNA space","Exclusively within tRNA genes","Expressed (# of datasets)","Maximum RPM")
tRFs_base<-tRFs[,1:3]
colnames(tRFs_base)<-c("Header","Sequence","Length")
tRFs_base$ID<-1:nrow(tRFs_base)
####for loop 
# tRFs_new<-data.frame(Header=NA,Sequence=NA,Length=NA)
# for (i in 1:nrow(tRFs_base)) {
#   temp=tRFs_base[grep(paste0('^',substr(tRFs_base[i,2],1,15)), tRFs_base$Sequence),]
#   if (nrow(temp)==1) {tRFs_new<-rbind(tRFs_new,temp)}
#   if (nrow(temp)>1) {
#     temp<-temp[order(-temp$Length),][1,]
#     tRFs_new<-rbind(tRFs_new,temp)
#   }
# }
# tRFs_new<-tRFs_new[!duplicated(tRFs_new$Header),]
# tRFs_new<-tRFs_new[!is.na(tRFs_new$Header),]
####foreach+dopar
registerDoParallel(cl<-makeCluster(20)) #Warning: check if your PC has 10 processors!
tRFs_new = foreach(i=1:nrow(tRFs_base),.combine=rbind)%dopar%{
  temp=tRFs_base[grep(paste0('^',substr(tRFs_base[i,2],1,15)), tRFs_base$Sequence),]
  if (nrow(temp)>1) {temp<-temp[order(-temp$Length),][1,]}
  temp
}
stopCluster(cl)
tRFs_new<-tRFs_new[,-3]
tRFs_new$ID<-1:nrow(tRFs_new)
tRFs_classified<-tRFs_new[!duplicated(tRFs_new$Header),]
writeFasta(tRFs_classified,"../Reference/tRFs/tRFs_classified.fa")
tRFs_merge<-merge(tRFs_base,tRFs_new,by="ID")
write.csv(tRFs_merge, "../Reference/tRFs_match.csv", row.names = F)

####ysRNA----
ysRNA_base<-readFasta("../Reference/ysRNA/ysRNA.fa")
ysRNA_base$Length<-nchar(ysRNA_base$Sequence)
ysRNA_base$ID<-1:nrow(ysRNA_base)
registerDoParallel(cl<-makeCluster(20)) #Warning: check if your PC has 10 processors!
ysRNA_new = foreach(i=1:nrow(ysRNA_base),.combine=rbind)%dopar%{
  temp=ysRNA_base[grep(paste0('^',substr(ysRNA_base[i,2],1,15)), ysRNA_base$Sequence),]
  if (nrow(temp)>1) {temp<-temp[order(-temp$Length),][1,]}
  temp
}
stopCluster(cl)
ysRNA_new<-ysRNA_new[,-3]
ysRNA_new$ID<-1:nrow(ysRNA_new)
ysRNA_classified<-ysRNA_new[!duplicated(ysRNA_new$Header),]
writeFasta(ysRNA_classified,"../Reference/ysRNA/ysRNA_classified.fa")
ysRNA_merge<-merge(ysRNA_base,ysRNA_new,by="ID")
write.csv(ysRNA_merge, "../Reference/ysRNA_match.csv", row.names = F)

####rsRNA----
rsRNA_base<-readFasta("../Reference/rsRNA/rsRNA.fa")
rsRNA_base$Length<-nchar(rsRNA_base$Sequence)
rsRNA_base$ID<-1:nrow(rsRNA_base)
registerDoParallel(cl<-makeCluster(20)) #Warning: check if your PC has 10 processors!
rsRNA_new = foreach(i=1:nrow(rsRNA_base),.combine=rbind)%dopar%{
  temp=rsRNA_base[grep(paste0('^',substr(rsRNA_base[i,2],1,15)), rsRNA_base$Sequence),]
  if (nrow(temp)>1) {temp<-temp[order(-temp$Length),][1,]}
  temp
}
stopCluster(cl)
rsRNA_new<-rsRNA_new[,-3]
rsRNA_new$ID<-1:nrow(rsRNA_new)
rsRNA_classified<-rsRNA_new[!duplicated(rsRNA_new$Header),]
writeFasta(rsRNA_classified,"../Reference/rsRNA/rsRNA_classified.fa")
rsRNA_merge<-merge(rsRNA_base,rsRNA_new,by="ID")
write.csv(rsRNA_merge, "../Reference/rsRNA_match.csv", row.names = F)
