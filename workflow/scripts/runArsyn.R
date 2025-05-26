log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(NOISeq)
library(readr)
library(dplyr)

##########################################
## ARSyN to reduce batch effect
##########################################

cat("Loading data\n")
load(snakemake@input[[1]])

mydata <- NOISeq::readData(
  data = full$M, 
  factors = full$target %>% select(group) %>% as.data.frame())

cat("Performing ARSyN for batch correction")
myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = FALSE)

##Saving everything
full <- list(M = assayData(myARSyN)$exprs, annot = full$annot, 
                              targets = full$targets)

cat("Saving RData \n") 
save(full, file=snakemake@output[[1]], compress="xz")