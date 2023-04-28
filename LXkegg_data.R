
setwd("D:/R-lin study/R packages/LXkegg")
library(openxlsx)

gene_data <- read.xlsx("gene_data.xlsx")
gene_data_log2FC <- read.xlsx("gene_data_log2FC.xlsx")

usethis::use_data(gene_data,overwrite = T)
usethis::use_data(gene_data_log2FC,overwrite = T)

rm(list=ls())

data(gene_data)
data(gene_data_log2FC)

