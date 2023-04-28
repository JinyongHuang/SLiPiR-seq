#Prepare a 'sample_barcode_list.txt' file in the working directory before start.
#Put the 'Reference' and 'Script' folder in previous directory.
mkdir tarfiles rawdata result_trimming result_mapping result_calling result_final
mkdir -p result_trimming/{index,cutadapter,splitbar,decompress,trimmomatic_15nt,trimmomatic_19nt,trimmomatic_23nt}
mkdir -p result_mapping/{rRNA,yRNA,miRNA,piRNA,hg38_15nt,hg38_23nt,Unmapped}
mkdir -p result_calling/{rsRNA,ysRNA,tRFs,miRNA,piRNA,mRNA,lncRNA,snoRNA,snRNA}

echo ----------------------------------------------------
echo Decompress and rename
mv *.tar tarfiles
for f in tarfiles/*.tar; do tar -xvf "$f"; done
mv *.fastq.gz rawdata
cd rawdata
ls *.fastq.gz > oldname.txt
ls *.fastq.gz | grep -o -E 'GDM[0-9]{8}-[A-Z]{1,}|R[1-2]'| paste -d "_" - - > newname.txt
awk 'NF{print $0 ".fastq.gz"}' newname.txt > newnames.txt
paste oldname.txt newnames.txt > rename.txt
awk -F '=>' '{print $1 $2}' rename.txt | mmv
cd ..
ls rawdata/*.fastq.gz | cut -d '/' -f 2 | cut -d '_' -f 1 | uniq > sample_list.txt
echo Complete

echo ----------------------------------------------------
echo Start geting index from R2
for i in `cat sample_list.txt`
do
zcat rawdata/${i}_R2.fastq.gz | awk 'NR%4==2{print(substr($0,1,8))}' > result_trimming/index/${i}_R2_index.txt
done
echo Complete

echo ----------------------------------------------------
echo Start cutting R1 adapter and polyA
date
cat sample_list.txt | parallel -j 80 cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o result_trimming/cutadapter/{}_cut_R1_temp.fastq.gz rawdata/{}_R1.fastq.gz
cat sample_list.txt | parallel -j 80 cutadapt -a AAAAAAAAAA -o result_trimming/cutadapter/{}_cut_R1.fastq.gz result_trimming/cutadapter/{}_cut_R1_temp.fastq.gz
rm result_trimming/cutadapter/*temp*
date
echo Complete

echo ----------------------------------------------------
echo Start spliting samples
date
cat sample_list.txt | parallel -j 12 python3 ../Script/barcode_extract.py result_trimming/cutadapter/{}_cut_R1.fastq.gz result_trimming/index/{}_R2_index.txt
#rm result_trimming/cutadapter/*.fastq.gz
#rm result_trimming/index/*.txt
date
echo Complete

echo ----------------------------------------------------
echo Start decompressing
date
for i in `cat sample_barcode_list.txt`
do
pigz -k -d result_trimming/splitbar/${i}_cut_R1.fastq.gz 
done
mv result_trimming/splitbar/*.fastq result_trimming/decompress
date
echo Complete

echo ----------------------------------------------------
echo Start trimming 15/19/23 nt
date
cat sample_barcode_list.txt | parallel -j 80 trimmomatic SE result_trimming/decompress/{}_cut_R1.fastq result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:15
cat sample_barcode_list.txt | parallel -j 80 trimmomatic SE result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq result_trimming/trimmomatic_19nt/{}_cut_trim_R1.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:19
cat sample_barcode_list.txt | parallel -j 80 trimmomatic SE result_trimming/trimmomatic_19nt/{}_cut_trim_R1.fastq result_trimming/trimmomatic_23nt/{}_cut_trim_R1.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:23
date
echo Complete