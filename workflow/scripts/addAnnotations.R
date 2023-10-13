log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
## Get the Work and Data dir
library(readr)
library(dplyr)
library(rtracklayer)

IS_XENA <- snakemake@params[["is_xena"]]
##############################################################################
## Get annotation data
## 1. Read the original annotation data
## 2. Read the new annotation data
## 3. Keep only genes/protein-coding present in original and new data
## 4. Query bioMart ensembl 80 for GC content 
## 5. Save data
###############################################################################
if(!IS_XENA) {
  cat("Getting annotation file \n")
  cat("Reading original file \n")
  annot <- rtracklayer::import(snakemake@input[["gdc_annot"]])
  annot <- as.data.frame(annot)
  
  ## Only protein coding genes
  annot <- annot %>% dplyr::select(gene_id, seqnames, start, end, width, type, 
                                   gene_type, gene_name) %>% 
    filter(type == "gene" & gene_type == "protein_coding")
  
  annot <- annot %>% mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1))) %>%
    select(-gene_type)

} else {
  annot <- read_tsv(snakemake@input[["xena_annot"]],
                    col_names  = c("gene_id", "gene_name", "seqnames", "start", "end", "strand"), 
                    skip = 1)
  annot <- annot %>% dplyr::mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)), 
                                   width = end - start)
  join_by <- "ensembl_id"
}

cat("Reading new file \n")
## Newest gencode file. April, 2021.
annot_new <-  rtracklayer::import(snakemake@params[["new_annot"]])
annot_new <- as.data.frame(annot_new)
annot_new <- annot_new %>% dplyr::select(gene_id, gene_name, type, gene_type) %>% 
  dplyr::filter(type == "gene" & gene_type == "protein_coding")

annot_new <- annot_new %>% dplyr::mutate(ensembl_id = unlist(lapply(strsplit(gene_id, "\\."), "[[", 1)))

## Keep genes that remain in the newest annotation file
## but get the newest names and keep only conventional chromosomes
## remove duplicates
annot <- annot %>% dplyr::select(-gene_name) %>%
  inner_join(annot_new %>% dplyr::select(ensembl_id, gene_name, gene_type), by = "ensembl_id") %>%
  filter(seqnames %in% paste0("chr", c(as.character(1:22), "X", "Y"))) %>% distinct()

cat('Annotation file new/old merge: ', paste(dim(annot), collapse=", "), '\n')

annot <- annot %>% mutate(chr = gsub("chr", "", seqnames)) %>%
  dplyr::rename(length = width) %>%
  dplyr::select(-seqnames) %>%   dplyr::select(gene_id, chr, everything())

cat('Annotation file. Final dimension: ', paste(dim(annot), collapse=", "), '\n')
save(annot, file=snakemake@output[["annot_rdata"]], compress="xz")
write_tsv(annot, snakemake@output[["annot_tsv"]])
cat('annot.RData saved \n')

##############################################################################
## Merging count and annotation
## 1. Build counts matrix. M=normal|tumour
## 2. Build targets matrix. targets=normal+tumor
## 3. Check M y targets integrity
## 4. Filter by annotation file 
## 5. Save the clean data
##############################################################################
{
  cat('Merging counts and annotations \n')
 
  normal_samples <- list(matrix=read_tsv(snakemake@input[["normal_matrix"]]),
                         targets=read_tsv(snakemake@input[["normal_targets"]])) 
  cancer_samples <- list(matrix=read_tsv(snakemake@input[["cancer_matrix"]]),
                         targets=read_tsv(snakemake@input[["cancer_targets"]])) 
  ## Raw counts
  M <- normal_samples$matrix %>% inner_join(cancer_samples$matrix, by = "gene_id")
  cat('Total number of features and samples: ', paste(dim(M), collapse=" ,"), '\n')
  
  ## Samples
  targets <- bind_rows(normal_samples$targets, cancer_samples$targets)
  
  ## Check M y targets integrity. Remove gene_ids col
  stopifnot(nrow(targets) == (ncol(M)-1))
  cat('Number of counts columns match sample number\n')
  
  ## Filter counts by annotation data
  cat('Adding biomart data\n')
  M <- M %>% semi_join(annot, by = "gene_id")
  annot <- annot %>% semi_join(M, by = "gene_id")
  
  cat('Total number of annotated (genes/protein-coding) features:', nrow(M), '\n')
  cat('Total number of samples:', ncol(M)-1, '\n')
  
  ## Save it as a matrix
  ids <- M$gene_id
  MM <- M %>% select(-gene_id) %>% as.matrix() 
  rownames(MM) <- ids
  MM <- MM[,targets$id]
  
  ## Make sure they are factors
  targets$group <- factor(targets$group, levels=c("cancer", "normal"))

  ##Save clean data
  cat('Saving raw full data \n')
  full <- list(M = MM, annot = annot, targets = targets)
  
  save(full, file=snakemake@output[["raw_rdata"]], compress="xz")
  cat("raw_full.RData saved \n")
}
