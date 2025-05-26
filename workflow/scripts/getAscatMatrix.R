log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(readr)
library(dplyr)
library(stringr)

tissue <- snakemake@params[["tissue"]]
cond <- snakemake@params[["type"]]
rawdir <- snakemake@params[["rawdir"]]

files_to_read <- list.files(path = rawdir, 
                            pattern = "\\.gene_level_copy_number.v36.tsv$", full.names = T, recursive = T)

cat("Got ", length(files_to_read), " samples\n")

if(length(files_to_read) == 0) {
  write_tsv(tibble(gene_id = c()), snakemake@output[[1]])
  write_tsv( tibble(id = c(), file = c(), file_id = c(), group = c()), snakemake@output[[2]])
  quit()
}
cat("Building matrices\n")
targets <- tibble(id = paste(tissue, cond, 1:length(files_to_read), sep = "_"), 
                      file = unlist(lapply(strsplit(files_to_read, "/"),  function(x) tail(x, n = 1))),
                      file_id = unlist(lapply(strsplit(files_to_read, "/"),  function(x) tail(x, n = 2)[1])),
                      group = cond)

all_data <- lapply(files_to_read, function(file) {
  id <- targets %>%
    filter(file == tail(strsplit(!!file, "/")[[1]], n = 1)) %>%
    pull(id)
  cnvs <- read_tsv(file)
  cnvs <- cnvs %>% 
    mutate(gene_id = as.character(lapply(str_split(gene_id, "\\."), "[[", 1))) %>%
    filter(!is.na(copy_number))
  return(cnvs %>% select(gene_id, !!id := copy_number))
})

all_data <- Reduce(inner_join, all_data) 
cat(nrow(all_data), " genes from ", ncol(all_data)-1,  " samples\n")

matrix <- all_data %>% 
  dplyr::select(gene_id, all_of(targets$id))

cat("Saving files\n")
write_tsv(matrix, snakemake@output[[1]])
write_tsv(targets, snakemake@output[[2]])
