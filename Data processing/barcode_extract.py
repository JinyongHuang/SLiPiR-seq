import re
import sys
import os
import gzip

R1=sys.argv[1]
R2=sys.argv[2]

R1_name = R1.split('/')[-1].split('_')[0].strip()
file_path=os.getcwd()

MB_R2=[]
R2 = open(R2,'r')
i=0
for line in R2:
	MB_R2.append(line.strip())
	i+=1
Total_reads = i
R2.close()

print ('Total reads in '+R1_name+' is: '+str(Total_reads))

R1=gzip.open(R1,'rb')
MB_uniq=['AGCGTAGC','TCGCCTTA','CTAGTACG','TTCTGCCT',
 	'GCTCAGGA','AGGAGTCC','CATGCCTA','GTAGAGAG',
 	'CCTCTCTG','CAGCCTCG','GTAAGGAG','TCCTCTAC',
 	'ACTGCATA','AAGGAGTA','CTAAGCCT','TCTCTCCG',
 	'GCGTAAGA','TAGATCGC','CTCTCTAT','AGAGTAGA']
for item in MB_uniq:
	exec ("R1_%s=gzip.open('./result_trimming/splitbar/'+str(R1_name+'_'+str(item)+'_cut_R1.fastq.gz'),'wt')"%item)
i=1
for line in R1:
	if i%4==1:
		R1_seq_ID=line.decode().split(' ')[0].split('/')[0].strip('\n')
	elif i%4==2:
		R1_seq=line.decode().strip('\n')
	elif i%4==0:
		R1_qua=line.decode().strip('\n')
		if(MB_R2[int(i/4-1)] == 'AGCGTAGC'):
			R1_AGCGTAGC.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'TCGCCTTA'):
			R1_TCGCCTTA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'CTAGTACG'):
			R1_CTAGTACG.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'TTCTGCCT'):
			R1_TTCTGCCT.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'GCTCAGGA'):
			R1_GCTCAGGA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'AGGAGTCC'):
			R1_AGGAGTCC.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'CATGCCTA'):
			R1_CATGCCTA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'GTAGAGAG'):
			R1_GTAGAGAG.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'CCTCTCTG'):
			R1_CCTCTCTG.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'CAGCCTCG'):
			R1_CAGCCTCG.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'GTAAGGAG'):
			R1_GTAAGGAG.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'TCCTCTAC'):
			R1_TCCTCTAC.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'ACTGCATA'):
			R1_ACTGCATA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'AAGGAGTA'):
			R1_AAGGAGTA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'CTAAGCCT'):
			R1_CTAAGCCT.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'TCTCTCCG'):
			R1_TCTCTCCG.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'GCGTAAGA'):
			R1_GCGTAAGA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'TAGATCGC'):
			R1_TAGATCGC.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'CTCTCTAT'):
			R1_CTCTCTAT.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
		elif(MB_R2[int(i/4-1)] == 'AGAGTAGA'):
			R1_AGAGTAGA.write(R1_seq_ID+'\n'+R1_seq+'\n'+'+'+'\n'+R1_qua+'\n')
	i+=1
	

			 











