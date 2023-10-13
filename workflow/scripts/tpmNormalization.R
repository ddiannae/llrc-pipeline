log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(readr)

cat("Reading raw.RData \n") 
load(snakemake@input[[1]])

geneLengths <- full$annot %>% 
  pull(length, name = "ensembl_id")

cat("Getting tpm Normalization")
rpk <- apply(full$M, 2, function(x) x/(geneLengths/1000))
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

tpm <- tpm[full$annot$gene_id, full$targets$id]

# Saving normalized counts
full <- list(M = tpm,
             annot = full$annot,
             targets = full$targets)

cat("Saving tpm_full.RData \n") 
save(full, file=snakemake@output[[1]], compress="xz")