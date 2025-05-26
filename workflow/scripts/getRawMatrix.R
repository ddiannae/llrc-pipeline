###############################################################################
##  LOSS OF LONG RANGE CO-EXPRESSION IN CANCER
##  Analysis of RNA-Seq data.
###############################################################################
## Diana Garcia - diana.gco@gmail.com
## Date: April, 2021
## 
## Original code. Dr. Cristobal Fresno - cristobalfresno@gmail.com
## Date:  2016-12-12
###############################################################################
## 1) Read Normal Data
## 2) Read Cancer Data
## 3) Read Annotation Data
## 4) Merge count and annotation
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
## Get the Work and Data dir
library(data.table)
library(readr)
library(dplyr)
library(janitor)
library(rtracklayer)

TISSUE <- snakemake@params[["tissue"]]
TYPE <-  snakemake@params[["type"]]
IS_XENA <- snakemake@params[["is_xena"]]
MCCORES <- as.numeric(snakemake@threads[[1]])
# ###############################################################################
# ## Get and check count matrices
# ## 1. Find files
# ## 2. Read data
# ## 4. Check sample sizes
# ## 5. Check genes order
# ## 6. Build target matrix
# ###############################################################################

if(!IS_XENA) {
  
  RAWDIR <- snakemake@params[["raw_dir"]]
 
  cat(paste0("Checking ", TYPE, " samples \n"))
     
  ## Find the files
  files_to_read <- list.files(path = RAWDIR, 
                              pattern = "\\.rna_seq.augmented_star_gene_counts.tsv$", full.names = T, recursive = T)
  
  cat("Got ", length(files_to_read), " samples\n")
  
  ## Build targets matrix
  targets <- tibble(id = paste(TISSUE, TYPE, 1:length(files_to_read), sep = "_"),
                    file = unlist(lapply(strsplit(files_to_read, "/"),  function(x) tail(x, n = 1))),
                    file_id = unlist(lapply(strsplit(files_to_read, "/"),  function(x) tail(x, n = 2)[1])),
                    group = TYPE)
  
  all_data <- lapply(files_to_read, function(file) {
    id <- targets %>%
      filter(file == tail(strsplit(!!file, "/")[[1]], n = 1)) %>%
      pull(id)
    data <- read_tsv(file, comment = "#") %>% 
      select(gene_id, gene_name, gene_type, !!id := unstranded) %>% 
      filter(!startsWith(gene_id, "N_"))
    
    return(data)
  })    
  
  ## Check samples sizes
  
  all_data <- Reduce(inner_join, all_data) 
  
  cat(nrow(all_data), " genes from ", ncol(all_data)-1,  " samples\n")
  
  all_data <- all_data %>%
    filter(gene_type %in% c("protein_coding")) %>%
    distinct(gene_id, .keep_all = T) %>%
    distinct(gene_name, .keep_all = T) 
  
 
  matrix <- all_data %>% 
    dplyr::select(gene_id, all_of(targets$id))
  
} else {
  XENA_COUNTS <- snakemake@input[[1]]
  XENA_SAMPLES <- snakemake@input[[2]]
  PRIMARY <- snakemake@params[["primary"]]
  
  matrix_samples <- read_tsv(XENA_SAMPLES) %>% clean_names()
  matrix_counts <-  fread(XENA_COUNTS, nThread = MCCORES)
  
  cat("Primary disease or tissue: ", PRIMARY, "\n")
  
  camelTissue <- paste0(toupper(substring(TISSUE, 1, 1)), substring(TISSUE, 2))
  
  cat("Tissue name: ", camelTissue, "\n")
  targets <- matrix_samples %>% dplyr::filter(primary_site == camelTissue &
                                                primary_disease_or_tissue == PRIMARY)
  cat("Got", nrow(targets), " samples\n")
  
  cat("Getting count matrix\n")
  ### Expected counts from XENA are in log2(expected_count+1)
  ### get them back to expected_count for downstream pipeline
  matrix <- matrix_counts %>% 
    dplyr::select_if(names(.) %in% c("sample", targets$sample)) %>%
    dplyr::select(sample, everything()) %>% 
    dplyr::mutate(across(-sample, ~ round(.x^2-1)))
  matrix[matrix < 0] <- 0
  
  cat(nrow(matrix), " genes from ", ncol(matrix)-1,  " samples\n")
    
  targets <- targets %>% dplyr::filter(sample %in% names(matrix)) %>%
    dplyr::mutate(id = paste(TISSUE, TYPE, 1:nrow(.), sep = "_"), group = TYPE) %>% 
    dplyr::rename(file = sample) %>% select(id, file, group)
    
  matrix <- matrix %>% 
    dplyr::select(gene_id = sample, all_of(targets$file))
  colnames(matrix) <- c("gene_id", targets$id)
}

cat(paste0("Matrices for ", TYPE, " ready.\n"))
cat("Saving matrix\n")
write_tsv(matrix, snakemake@output[["matrix"]])
cat("Saving samples\n")
write_tsv(targets, snakemake@output[["samples"]])
