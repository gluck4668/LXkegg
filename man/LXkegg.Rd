\name{LXkegg}
\alias{LXkegg}
\title{KEGG analysis}
\usage{ LXkegg(gene_file,group1,group2,species)}
\description{
KEGG pathways and GO enrichment analyses}
\examples{

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

setwd("D:/R_example/LXkegg_example")
gene_file="gene_data_log2FC.xlsx"
group1="Normal"
group2="Model"
species= "mouse"  # The species should be "human", "mouse", or "rat"

LXkegg(gene_file,group1,group2,species)

}
