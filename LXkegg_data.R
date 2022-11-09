
setwd("D:/R-lin study/R packages/LXkegg")
library(openxlsx)

gene_data_example <- read.xlsx("gene_data_example.xlsx")

usethis::use_data(gene_data_example,overwrite = T)

rm(list=ls())

data(gene_data_example)

