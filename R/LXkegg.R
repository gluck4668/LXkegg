
#----------------KEGG analysis--------------------------------------#

LXkegg <- function(gene_file,group1,group2,species){


  kegg_df <- data_processing(gene_file,group1,group2,species)

  kegg_pathways <- kegg_pathways(kegg_df)

  GO_enrich <- GO_enrich(kegg_df)

  kegg_pathways <- kegg_pathways$kegg_pathways

  GO_enrich <- GO_enrich$GO_enrich

  dir.file <- kegg_df$dir.file

  p <- GO_enrich+kegg_pathways+plot_layout(ncol=2,nrow=1,widths = c(3,1))

  GO_kegg <- paste(species,"GO enrichment and kegg pathways.png")
  GO_kegg <-paste0(dir.file,"/",GO_kegg)
  ggsave(GO_kegg, p, width=2400, height =1000, dpi=150,units = "px")

  dev.new(width=18,height=7, noRStudioGD = T)

  print("The analysis results can be found in the folder of <analysis result>")

  p

}
