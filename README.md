## Introduction
The codes within the ‘SLiPiR-seq’ repository are used for data processing and data analysis of sequencing data produced by the SLiPiR-seq cell-free RNA profiling technology. The codes in the ‘Data processing’ repository start with raw sequencing data in fastq format from Ilumina platform and end with the calling of read count for nine different RNA categories. The code in the ‘Data analysis’ repository are used for quality control of each sample, differential expression analysis between different groups of samples, selection of discriminating features, and establishment of machine learning models.
## Requirements
Computer running a Linux system Cluster computing is highly recommended when working with fastq/bam files. All computational analyses were performed using Linux shell, Python 3, or R 4.2.\
Linux shell software: cutadapt (version 2.10); trimmomatic (version 0.39); bowtie2 (version 2.4.4); subread package (version 2.0.3); \
R packages:  tidyverse (version 2.0.0); ggpubr (version 0.6.0); ggplot2 (version 3.4.2); ggpmisc (version  0.5.2); ggrepel (version 0.9.3); RColorBrewer (version 1.1-3); patchwork (version 1.1.2); DESeq2 (version 1.34.0); EnhancedVolcano (version 1.12.0); caret (version 6.0-94); foreach (version 1.5.2); doParallel (version 1.0.17); glmnet (version  4.1-7); randomForest (version 4.7-1.1); ROCR (version 1.0-11);  pROC (version 1.18.0); reshape2 (version 1.4.4)
## Procedure
1.	Raw data FASTQ files to processed FASTQ files (clean reads).\
```./Data processing/1. Cut_splitbarcode_decompress.sh```
2.	(optional) Generating reference genome for rsRNA, tsRNA and ysRNA. Or skip this step and use the provided reference genomes directly.\
```./Data processing/2. rsysRNA_db_update.sh```
3.	Read calling for nine RNA categories.\
```./Data processing/3. RNA read count calling.sh```
4.	Quality control and normalization.\
```Rscript ./Data analysis/1.Quality control and normalization/1.Raw read counts.R ```\
```Rscript ./Data analysis/1.Quality control and normalization/2.Normalization-RPM-colsum.R ```
5.	Case and control study – One cancer type\
Differential expression analysis:\
```Rscript ./Data analysis/2.Lung cancer/1.DE/1.DE_LC_SZDE_vs_NOR_SZBA.R```\
```Rscript ./Data analysis/2.Lung cancer/1.DE/2.Candidates_reads_boxplot.R```\
```Rscript ./Data analysis/2.Lung cancer/1.DE/2.Heatmap.R```\
Feature selection:\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/1.Feature selection-Boruta.R```\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/1.Feature selection-LASSO.R```\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/1.Feature selection-TopN.R```\
Machine learning analysis:\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/2.AUC single RNA.R```\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/2.ROC single.R```\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/2.goi_reads_boxplot.R```\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/3.Combination test.R```\
```Rscript ./Data analysis/2.Lung cancer/2.Machine learning/4.Best combination.R```
6.	Case and control study – Multiple cancer types.\
Differential expression analysis:\
```Rscript ./Data analysis/3.Multiple cancer types/1.DE/DE_BRC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/1.DE/DE_CRC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/1.DE/DE_GC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/1.DE/DE_HCC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/1.DE/DE_LC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/1.DE/DE_cancers_vs._cancerfree.R```\
Machine learning analysis:\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/ML_BRC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/ML_CRC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/ML_GC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/ML_HCC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/ML_LC_vs._others.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/ML_cancers_vs._cancerfree.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/Heatmap.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/tsne.R```\
```Rscript ./Data analysis/3.Multiple cancer types/2.Machine learning/Validation cohort.R```\
