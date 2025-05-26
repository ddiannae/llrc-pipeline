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

tpm_cancer <- read_tsv(snakemake@input[["cancer"]])
tpm_normal <- read_tsv(snakemake@input[["normal"]])

tpm <- tpm_cancer %>%
  full_join(tpm_normal, by = "gene")

genes <- tpm$gene
tpm <- as.matrix(tpm %>% select(-gene))
row.names(tpm) <- genes

print("Running cibersort")

bs <- as.data.frame(CIBERSORT(sigmatrix, tpm, perm = 100, QN = FALSE)) %>%
  clean_names() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(snakemake@output[["cibersort"]])

print("Running xcell")

res_xcell = as.data.frame(deconvolute(tpm, "xcell"))
rownames(res_xcell) <- res_xcell$cell_type

res_xcell %>%
  select(-cell_type) %>% 
  t() %>%
  clean_names() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(snakemake@output[["xcell"]])
