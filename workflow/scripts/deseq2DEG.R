log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(DESeq2)
library(dplyr)
library(ggplot2)
library(readr)

cat("Reading raw.RData \n") 
load(snakemake@input[["raw"]])

raw_full <- full

cat("Reading deseq2.RData \n") 
load(snakemake@input[["deseq2"]])

## The deseq2.Rdata has vst transformed counts. They shouldn't 
## be used for DEG but we do wanna extract the samples that 
## remained after the outlier removal 

full$targets %>% nrow()

targets <- raw_full$targets %>% 
  filter(id %in% full$targets$id)
M <- raw_full$M[ ,targets$id]

genes <- raw_full$annot %>%
  pull(ensembl_id, name = gene_id)
M <- M[names(genes), ]
rownames(M) <- genes

targets$group <- relevel(targets$group, ref = "normal")

dds <- DESeqDataSetFromMatrix(countData = M,
                              colData = targets,
                            design = ~group)

dds <- DESeq(dds)

resLFC <- lfcShrink(dds, coef="group_cancer_vs_normal", type="apeglm")

as_tibble(resLFC, rownames = "ensembl_id") %>%
  mutate(padj = ifelse(padj > 0.05, 1, padj), 
         log2FC = ifelse(padj > 0.05, 0, log2FoldChange)) %>%
  select(-log2FoldChange) %>%
  write_tsv(snakemake@output[["deg_results"]])
