log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(readr)
library(immunedeconv)
library(janitor)
source(snakemake@params[["cibersort"]])

outdir <- snakemake@params[["outdir"]]

if (!dir.exists(outdir)) {
  dir.create(outdir)
}

sigmatrix <- snakemake@params[["signature_matrix"]]

tpm <- read_tsv(snakemake@input[[1]])
genes <- tpm$gene
tpm <- as.matrix(tpm %>% select(-gene))
row.names(tpm) <- genes

print("Running cibersort")

bs <- as.data.frame(CIBERSORT(sigmatrix, tpm, perm = 100, QN = FALSE)) %>%
  clean_names() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(snakemake@output[[1]])
