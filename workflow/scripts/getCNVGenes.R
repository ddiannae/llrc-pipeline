#######################
## This script extracts the mean copy number values for each gene
##  from an ASCAT matrix.
#######################

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
     
library(dplyr)
library(tibble)
library(vroom)

cat("Reading matrix\n")
ascat_matrix <- vroom(snakemake@input[["ascat_matrix"]])

if(nrow(ascat_matrix > 1)) {
  ascat_matrix <- ascat_matrix %>% 
    column_to_rownames("gene_id") %>% 
    as.matrix()
  
  cat("Getting CNVs\n")
  cnv_means <- rowMeans(ascat_matrix)
  cnv_means <- tibble(gene_id = names(cnv_means), mean_val = cnv_means)
  cat("Writing matrix\n")
} else {
  cnv_means <- tibble()
}

vroom_write(cnv_means, file = snakemake@output[["cnv_genes"]])

