log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(tibble)
library(dplyr)
library(DESeq2)

cat("Reading raw.RData \n") 
load(snakemake@input[[1]])

cat("Getting deseq2 Normalization")
# Normalization with DESeq2
dds <- DESeqDataSetFromMatrix(countData = full$M,
                              colData = column_to_rownames(full$targets, var = "id"), 
                              design = ~group)
dds <- estimateSizeFactors(dds)

vst <- vst(dds, blind=TRUE)
vst_counts <- assay(vst)
vst_counts <- vst_counts[full$annot$gene_id, full$targets$id]

# Saving normalized counts
full <- list(M = vst_counts,
     annot = full$annot,
     targets = full$targets)

cat("Saving deseq2_full.RData \n") 
save(full, file=snakemake@output[[1]], compress="xz")
saveRDS(dds, snakemake@output[["dds"]])
