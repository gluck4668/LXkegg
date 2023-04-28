\name{LXkegg}
\alias{LXkegg}
\title{KEGG analysis}
\usage{ LXkegg(gene_file,group1,group2,species)}
\description{
KEGG pathways and GO enrichment analyses}
\examples{

install.packages("devtools")
library(devtools)

install_github("gluck4668/LXkegg",force=T)

library(LXkegg)

??LXkegg
#---------------------
data(gene_data)
data(gene_data_log2FC)

#--------------------

rm(list=ls())

setwd("D:/桌面/R_example/LXkegg_example")

gene_file="gene_data_log2FC.xlsx" # The data can be a list of gene SYMBOL, or two list including gene SYMBOL and log2FoldChange

group1="Didymin"
group2="Model"

species= "mouse"  # The species should be "human", "mouse", or "rat"

LXkegg(gene_file,group1,group2,species)

}
