log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

cat("Reading files\n")

cnv_matrix <- vroom(snakemake@input[[1]])
if(nrow(cnv_matrix) == 0) {
  vroom_write(tibble(),snakemake@output[["pcmatrix"]])
  vroom_write(tibble(),snakemake@output[["heatmap"]])
  quit()
}
genes <- unlist(cnv_matrix[, 1], use.names = F)
cnv_matrix <- as.matrix(cnv_matrix[,-1])
rownames(cnv_matrix) <- genes

annot <- vroom(snakemake@params[["biomart"]], 
               col_names = c("id", "chr", "start", "end", "gc_content", "type"),
               skip = 1)
chrs <- c(as.character(1:22), "X", "Y")

cat("Merging annotations\n")
annot <- annot %>% 
  filter(id %in% genes, chr %in% chrs) %>%  
  mutate(chr = factor(chr, levels = chrs)) %>%  
  arrange(chr, start) %>%
  filter(type == "protein_coding")
 
genes <- annot$id
cnv_matrix <- cnv_matrix[genes, ]

as_tibble(cnv_matrix, rownames = "gene_id") %>%
  select(gene_id, everything()) %>%
  vroom_write(snakemake@output[["pcmatrix"]])

cat("Building heatmap\n")
col_fun = colorRamp2(c(0, 2, 6), c("blue", "white", "red"))
chromosomes.pal <- c("#D909D1", "#0492EE", "#5DA0CB", "#106F35", "#5BD2AE", "#199F41", 
                     "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
                     "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#800092", "#B925AE",
                      "#B1B719", "#CB97E8", "#130B9E", "#E12B29", "#79A5B9", "#C8CDCE")

names(chromosomes.pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                            "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21", "Y")

chrs <- as.data.frame(annot %>% select(chr) %>% 
                        rename(Chr = chr))

ht_opt("raster_temp_image_max_width" = 15000, 
       "raster_temp_image_max_height" = 15000) 
       

ha <- HeatmapAnnotation(df = chrs, name = "Chr", show_annotation_name = F,
                        col = list(Chr = chromosomes.pal),
                        which = "row", width = unit(0.5, "cm"))

ht <- Heatmap(cnv_matrix, cluster_rows = F, cluster_columns = T, 
        col = col_fun, name = "CN", na_col = "#000000",
        show_column_names = F, show_row_names = F, use_raster = TRUE,
        raster_by_magick = TRUE)

cat("Saving heatmap\n")
png(filename = snakemake@output[["heatmap"]], width = 600, height = 1200)
ht + ha
dev.off()

ht_opt(RESET = TRUE)
