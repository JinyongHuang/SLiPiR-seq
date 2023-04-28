import os
import re
import sys

R1 = sys.argv[1]
outputname = sys.argv[2]
file_path=os.getcwd()
d_ysR={}
d_cal={}
ysR=open('../Reference/ysRNA/ysRNA.fa','r')
ysR_cal_output = open(file_path+'/result_calling/ysRNA/'+str(outputname)+'ysRNA.txt','a')
i=1
for line in ysR:
        if i%2==1:
                ID = line.split('\t')[0].strip().split('>')[1].strip()
        if i%2==0:
                seq = line.strip()
                d_ysR[seq]=ID
                d_cal[ID]=0
        i+=1
R1=open(R1,'r')
i=1
for line in R1:
        if i%4==1:
                R1_seq_ID=line.split(' ')[0].strip('\n')
        elif i%4==2:
                R1_seq=line.strip('\n')
        elif i%4==0:
                R1_qua=line.strip('\n')
                if R1_seq in d_ysR.keys():
                        ID = d_ysR[R1_seq].strip()
                        d_cal[ID]+=1
        i+=1
for items in d_cal.keys():
        ysR_cal_output.write(items+'\t'+str(d_cal[items])+'\n')
