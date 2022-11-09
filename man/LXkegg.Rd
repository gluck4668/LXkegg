\name{LXkegg}
\alias{LXkegg}
\title{KEGG analysis}
\usage{ LXkegg(gene_file,group1,group2)}
\description{
KEGG pathways and GO enrichment analyses}
\examples{
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXkegg")

library(LXkegg)

#---------------------
data(gene_data_example)

#--------------------

rm(list=ls())

setwd("D:/R-lin study/R packages/LXkegg/LXkegg_example")
gene_file="gene_data.xlsx"
group1="Normal"
group2="Model"

LXkegg(gene_file,group1,group2)

}
