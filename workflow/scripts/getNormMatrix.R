#######################################
### Get normal and cancer matrix from norm data
#######################################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(readr)

cat("Reading norm data \n") 
load(snakemake@input[[1]])

id <- snakemake@params[["gene_id"]]

if(id == "name") {
  genes <- full$annot %>%
    pull(gene_name, name = gene_id)
} else {
  genes <- full$annot %>%
    pull(ensembl_id, name = gene_id)
}

genes <- genes[!duplicated(genes)]
M <- full$M[names(genes), ]
rownames(M) <- genes

M[, full$targets %>%
    filter(group == "normal") %>%
    pull(id)] %>%
  as_tibble(rownames = "gene") %>%
  select(gene, everything()) %>%
  write_tsv(snakemake@output[["normal"]])

M[, full$targets %>%
    filter(group == "cancer") %>%
    pull(id)] %>%
  as_tibble(rownames = "gene") %>%
  select(gene, everything()) %>%
  write_tsv(snakemake@output[["cancer"]])

tibble(genes) %>%
  write_tsv(snakemake@output[["genes"]], col_names = FALSE)
