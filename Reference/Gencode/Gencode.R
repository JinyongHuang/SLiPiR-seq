library(rtracklayer)
library(plyranges)

df<-import("gencode.v41.annotation.gtf")
unique(df$gene_type)

protein_coding<- df %>% filter(gene_type == "protein_coding" | 
                                 gene_type == "IG_C_gene" |
                                 gene_type == "IG_D_gene" |
                                 gene_type == "IG_J_gene" |
                                 gene_type == "IG_V_gene" |
                                 gene_type == "TR_C_gene" |
                                 gene_type == "TR_D_gene" |
                                 gene_type == "TR_J_gene" |
                                 gene_type == "TR_V_gene")
export(protein_coding,"gencode.v41.protein_coding.gtf")

lncRNA<-df %>% filter(gene_type == "lncRNA")
export(lncRNA,"gencode.v41.lncRNA.gtf")

snRNA<-df %>% filter(gene_type == "snRNA")
export(snRNA,"gencode.v41.snRNA.gtf")

snoRNA<-df %>% filter(gene_type == "snoRNA")
export(snoRNA,"gencode.v41.snoRNA.gtf")
