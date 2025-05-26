log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)
library(ComplexHeatmap)
library(circlize)


t <- snakemake@params[["tissue"]]
color_pal <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")
names(color_pal) <- labels

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

ciber <- vroom(snakemake@input[["cibersort"]]) %>%
  mutate(type = if_else(grepl("cancer", sample), "Cancer", "Normal"))
  
mciber <- ciber %>%
  select(cd10, cd31, cd45, epcam) %>%
  as.matrix()

rownames(mciber) <- ciber %>%
  pull(sample)

conds <- ciber %>% pull(type, name = sample)
col_fun = colorRamp2(c(0, 1), c("white", "navy"))

ht <- Heatmap(mciber, cluster_rows = TRUE, cluster_columns =  FALSE, show_row_names = F, 
        show_column_names = TRUE, name = "Score", col = col_fun,
        column_labels = c("CD10", "CD31", "CD45", "EPCAM"), 
        column_names_rot = 0,
        column_names_gp =  gpar(fontsize = 20),
        column_names_centered = T,
        left_annotation = rowAnnotation(Condition = conds, 
                                        col = list(Condition = color_pal), 
                                        show_annotation_name = F, 
                                        annotation_legend_param = list(title_gp = gpar(fontsize = 20), 
                                        labels_gp = gpar(fontsize = 18))), 
        heatmap_legend_param = list(legend_height = unit(3, "cm"), 
                                    title_gp = gpar(fontsize = 20), 
                                    labels_gp = gpar(fontsize = 18), 
                                    title_position = "leftcenter-rot"),
        column_title = paste0(capitalize(t), " Cibersort"),
        column_title_gp = grid::gpar(fontsize = 25)) 

png(snakemake@output[["cibersort"]], width = 500, height = 1000)
draw(ht,   padding = unit(0.5, "cm"))
dev.off()


xcell <- vroom(snakemake@input[["xcell"]])
xcell <- xcell %>%
mutate(type = if_else(grepl("cancer", sample), "Cancer", "Normal"))  %>%
  select(-immune_score, -microenvironment_score,  -stroma_score)


mxcell <- xcell %>%
  select(-sample, -type) %>%
  as.matrix()

rownames(mxcell) <- xcell %>%
  pull(sample)

conds <- xcell %>% pull(type, name = sample)
col_fun = colorRamp2(c(0, max(mxcell)), c("white", "navy"))

ht <- Heatmap(mxcell, cluster_rows = TRUE, cluster_columns =  FALSE, show_row_names = F, 
              show_column_names = TRUE, name = "Score", col = col_fun,
              column_names_rot = 90,
              column_names_centered = F,
              left_annotation = rowAnnotation(Condition = conds, 
                                              col = list(Condition = color_pal), 
                                              show_annotation_name = F, 
                                              annotation_legend_param = list(title_gp = gpar(fontsize = 20), 
                                                                             labels_gp = gpar(fontsize = 18))), 
              heatmap_legend_param = list(legend_height = unit(3, "cm"), 
                                          title_gp = gpar(fontsize = 20), 
                                          labels_gp = gpar(fontsize = 18), 
                                          title_position = "leftcenter-rot"),
              column_title = paste0(capitalize(t), " xCell"),
              column_title_gp = grid::gpar(fontsize = 25)) 

 png(snakemake@output[["xcell"]], width = 1000, height = 500)
 draw(ht,   padding = unit(0.5, "cm"))
 dev.off()
