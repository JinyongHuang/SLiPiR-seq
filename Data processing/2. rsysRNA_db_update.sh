echo ----------------------------------------------------
echo Start updating rsRNA database
date
cat sample_barcode_list.txt | parallel -j 80 bowtie2 --score-min C,0 --no-unal -x ../Reference/rRNA/Homo_rRNA -S result_mapping/rRNA/{}_cut_trim_map2rRNA.sam -U result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq '&>' result_mapping/rRNA/{}_cut_trim_map2rRNA.log
doawk() {
  INPUT=$1
  awk '{if ($6~/^[0-9]*M$/) print $3 "\t" $10}' "$INPUT" | sort | uniq -c | sort -nr
} 
export -f doawk
cat sample_barcode_list.txt | parallel -j 80 doawk result_mapping/rRNA/{}_cut_trim_map2rRNA.sam '>' result_mapping/rRNA/{}_tmp.txt
doawk() {
  INPUT=$1
  awk 'BEGIN{i=1}{if($1>10 && $3~/^[ATCG]*$/){print ">"$2"-"i++"\n"$3;}}' "$INPUT" 
} 
export -f doawk
cat sample_barcode_list.txt | parallel -j 80 doawk result_mapping/rRNA/{}_tmp.txt '>' result_mapping/rRNA/{}_tmp.fa
for i in `cat sample_barcode_list.txt`
do
python3 ../Script/rsRNA_db_update.py result_mapping/rRNA/${i}_tmp.fa
done
rm result_mapping/rRNA/*tmp*
#rm result_mapping/rRNA/*.sam
date
echo Complete

echo ----------------------------------------------------
echo Start updating ysRNA database
date
cat sample_barcode_list.txt | parallel -j 80 bowtie2 --score-min C,0 --no-unal -x ../Reference/yRNA/Homo_yRNA -S result_mapping/yRNA/{}_cut_trim_map2yRNA.sam -U result_trimming/trimmomatic_15nt/{}_cut_trim_R1.fastq '&>' result_mapping/yRNA/{}_cut_trim_map2yRNA.log
doawk() {
  INPUT=$1
  awk '{if ($6~/^[0-9]*M$/) print $3 "\t" $10}' "$INPUT" | sort | uniq -c | sort -nr
} 
export -f doawk
cat sample_barcode_list.txt | parallel -j 80 doawk result_mapping/yRNA/{}_cut_trim_map2yRNA.sam '>' result_mapping/yRNA/{}_tmp.txt
doawk() {
  INPUT=$1
  awk 'BEGIN{i=1}{if($1>10 && $3~/^[ATCG]*$/){print ">"$2"-"i++"\n"$3;}}' "$INPUT" 
} 
export -f doawk
cat sample_barcode_list.txt | parallel -j 80 doawk result_mapping/yRNA/{}_tmp.txt '>' result_mapping/yRNA/{}_tmp.fa
for i in `cat sample_barcode_list.txt`
do
python3 ../Script/ysRNA_db_update.py result_mapping/yRNA/${i}_tmp.fa
done
rm result_mapping/yRNA/*tmp*
#rm result_mapping/yRNA/*.sam
date
echo Complete

echo ----------------------------------------------------
echo Start classifying tRFs, rsRNA and ysRNA references
date
Rscript ../Script/Reference.R
date
echo Complete
