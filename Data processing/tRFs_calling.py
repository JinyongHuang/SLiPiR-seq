import os
import re
import sys

R1 = sys.argv[1]
outputname = sys.argv[2]
file_path=os.getcwd()
d_ID={}
d_tRFID={}
tRFs=open('../Reference/tRFs/tRFs_original.fa','r')
for line in tRFs:
	ID = line.split('\t')[0].strip()
	seq = line.split('\t')[1].strip()
	d_ID[seq]=ID
for key in d_ID.keys():
	d_tRFID[d_ID[key]]=0
R1_name_tRF_ID = str(R1.split('/')[-1].split('.')[0].split('_R1')[0]+'_tRFs.txt')
R1_output_tRF_ID=open(file_path+'/result_calling/tRFs/'+R1_name_tRF_ID,'a')
R1=open(R1,'r')
i=1
for line in R1:
	if i%4==1:
		R1_seq_ID=line.split(' ')[0].strip('\n')
	elif i%4==2:
		R1_seq=line.strip('\n')
	elif i%4==0:
		R1_qua=line.strip('\n')
		if R1_seq in d_ID.keys():
			d_tRFID[d_ID[R1_seq]]+=1
	i+=1
for ID in d_tRFID.keys():
	R1_output_tRF_ID.write(str(ID)+'\t'+str(d_tRFID[ID])+'\n')

