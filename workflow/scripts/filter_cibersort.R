log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(dplyr)

ADJ <-  snakemake@input[[1]]
matrix <- vroom(ADJ)

matrix$nas <- rowSums(is.na(matrix))

matrix <- matrix %>%
  filter(nas <= ncol(matrix) /2) %>%
  select(-nas)

matrix$means <- rowMeans(matrix %>%
                           select(where(is.numeric)))
matrix <- matrix %>% 
  filter(means >= 10) %>%
  select(- means)
  
vroom_write(matrix, snakemake@output[[1]])