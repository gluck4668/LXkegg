
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXkegg")

library(LXkegg)

??LXkegg
#---------------------
data(gene_data)
data(gene_data_log2FC)

#--------------------

rm(list=ls())

#devtools::load_all()

gene_file="gene_data.xlsx"
gene_file="gene_data_log2FC.xlsx" # The data can be a list of gene SYMBOL, or two list including gene SYMBOL and log2FoldChange


group1="Treatment"
group2="Model"

species= "rat"  # The species should be "human", "mouse", or "rat"

LXkegg(gene_file,group1,group2,species)




