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

echo ----------------------------------------------------
echo Start calling rsRNA
date
cat sample_barcode_list.txt | parallel -j 80 python3 ../Script/rsRNA_calling.py result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq {}_
date
echo Complete

echo ----------------------------------------------------
echo Start calling ysRNA
date
cat sample_barcode_list.txt | parallel -j 80 python3 ../Script/ysRNA_calling.py result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq {}_
date
echo Complete

echo ----------------------------------------------------
echo Start calling tRFs
date
cat sample_barcode_list.txt | parallel -j 80 python3 ../Script/tRFs_calling.py result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq {}_
date
echo Complete

echo ----------------------------------------------------
echo Start calling piRNAs
date
cat sample_barcode_list.txt | parallel -j 80 bowtie2 -x ../Reference/piRNA/piR_human -S result_mapping/piRNA/{}_map2piRNA.sam -U result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq '&>' result_mapping/piRNA/{}_map2piRNA.log
doawk() {
  INPUT=$1
  awk '{if ($6~/^[0-9]*M$/) print $3}' "$INPUT" | sort | uniq -c | sort -nr
} 
export -f doawk
cat sample_barcode_list.txt | parallel -j 80 doawk result_mapping/piRNA/{}_map2piRNA.sam '>' result_calling/piRNA/{}_piRNA.txt
rm result_mapping/piRNA/*.sam
date
echo Complete

echo ----------------------------------------------------
echo Start calling miRNA
date
cat sample_barcode_list.txt | parallel -j 80 bowtie2 -x ../Reference/miRNA/hsa_miRNAs -S result_mapping/miRNA/{}_map2miRNA.sam -U result_trimming/trimmomatic_19nt/{}_cut_trim_R1.fastq '&>' result_mapping/miRNA/{}_map2miRNA.log
doawk() {
  INPUT=$1
  awk '{if ($6~/^[0-9]*M$/) print $3}' "$INPUT" | sort | uniq -c | sort -nr
} 
export -f doawk
cat sample_barcode_list.txt | parallel -j 80 doawk result_mapping/miRNA/{}_map2miRNA.sam '>' result_calling/miRNA/{}_miRNA.txt
rm result_mapping/miRNA/*.sam
date
echo Complete

echo ----------------------------------------------------
echo Start mapping to hg38
##3.25GB (1.3%) memory per sample
date
cat sample_barcode_list.txt | parallel -j 60 bowtie2 -x ../Reference/GRCh38/hg38 -U result_trimming/trimmomatic_23nt/{}_cut_trim_R1.fastq -S result_mapping/hg38_23nt/{}_map2hg38.sam '&>' result_mapping/hg38_23nt/{}_map2hg38.log
cat sample_barcode_list.txt | parallel -j 60 bowtie2 -x ../Reference/GRCh38/hg38 -U result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq -S result_mapping/hg38_15nt/{}_map2hg38.sam '&>' result_mapping/hg38_15nt/{}_map2hg38.log
date
echo Complete

echo ----------------------------------------------------
echo Start calling mRNA, lncRNA, snoRNA, and snRNA
date
featureCounts -T 60 -O -s 1 -t transcript -g gene_name -a ../Reference/Gencode/gencode.v41.protein_coding.gtf -o result_calling/mRNA/mRNA_S1_transcript.txt result_mapping/hg38_23nt/*.sam
featureCounts -T 60 -O -s 1 -t transcript -g gene_name -a ../Reference/Gencode/gencode.v41.lncRNA.gtf -o result_calling/lncRNA/lncRNA_S1_transcript.txt result_mapping/hg38_23nt/*.sam
featureCounts -T 60 -O -s 1 -t transcript -g gene_name -a ../Reference/Gencode/gencode.v41.snRNA.gtf -o result_calling/snRNA/snRNA_S1_transcript.txt result_mapping/hg38_15nt/*.sam
featureCounts -T 60 -O -s 1 -t transcript -g gene_name -a ../Reference/Gencode/gencode.v41.snoRNA.gtf -o result_calling/snoRNA/snoRNA_S1_transcript.txt result_mapping/hg38_15nt/*.sam
rm result_mapping/hg38/*23nt.sam
date
echo Complete

echo ----------------------------------------------------
echo Start calculating total/trimmed/unmapped reads
date
find . -name "*_reads.txt" -delete
cat sample_barcode_list.txt | parallel -j 80 samtools view -f 4 result_mapping/hg38_15nt/{}_map2hg38.sam '>' result_mapping/Unmapped/{}_unmapped.sam
cat sample_barcode_list.txt | parallel -j 80 wc -l result_mapping/Unmapped/{}_unmapped.sam '>>' result_final/Unmapped_reads.txt
cat sample_barcode_list.txt | parallel -j 80 wc -l result_trimming/decompress/{}_cut_R1.fastq | awk '{print $1/4"\t"$2}' >> result_final/Total_reads.txt
cat sample_barcode_list.txt | parallel -j 80 wc -l result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq | awk '{print $1/4"\t"$2}' >> result_final/Trimmed_reads.txt
date
echo Complete

echo ----------------------------------------------------
echo Start making counts matrix
date
Rscript ../Script/Counts_matrix.R
date
echo Complete

tar cf result_final.tar result_final/*
#sz result_final.tar 

# echo ----------------------------------------------------
# echo Optional: sam to bam
# date
# cat sample_barcode_list.txt | parallel -j 40 samtools view -S result_mapping/hg38_15nt/{}_map2hg38.sam -b '>' result_mapping/hg38_15nt/{}_map2hg38.bam
# rm result_mapping/hg38_15nt/*.sam
# cat sample_barcode_list.txt | parallel -j 40 samtools view -S result_mapping/Unmapped/{}_unmapped.sam -b '>' result_mapping/Unmapped/{}_unmapped.bam
# rm result_mapping/Unmapped/*.sam
# date
# echo Complete