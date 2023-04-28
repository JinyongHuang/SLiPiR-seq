import re
import sys
import os

tmp_db=sys.argv[1]
file_path=os.getcwd()
db=open('../Reference/ysRNA/ysRNA.fa','r')
output=open('../Reference/ysRNA/ysRNA.fa','a')
tmp_db=open(tmp_db,'r')
d={}
i=0
for line in db:
	i+=1
	if i%2==1:
		ID = line.strip()
	if i%2==0:
		seq = line.strip()
		d[seq]=ID
i=int(i/2)
y=0
for line in tmp_db:
	y+=1
	if y%2==1:
		inf = str(line.split('-')[0]+'-')
	if y%2==0:
		seq = line.strip()
		if seq not in d.keys():
			i=i+1
			output.write(str(inf+str(i))+'\n'+seq+'\n')

		 
