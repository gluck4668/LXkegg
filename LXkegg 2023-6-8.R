
#----------------KEGG analysis--------------------------------------#

LXkegg <- function(gene_file,group1,group2,species){

  #list all the packages that have been installed
  all_packages <- data.frame(installed.packages())

  #To judge whether a package was installed. If not, it will be installed.
  pack <- data.frame(c("devtools","BiocManager","ggnewscale","R.utils",
                       "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych","ggplot2",
                       "ggrepel","RColorBrewer", "ggthemes","rticles","grid","patchwork","Hmisc","pak") )

  # To judge whether a package was included in the all_packages: %in%
  pack$type <- pack[,1] %in% all_packages$Package

  for (i in 1:nrow(pack)){
    if(pack[i,2]==FALSE)
      install.packages(pack[i,1],update = F,ask = F)
  }
  rm(i)

  # 批量library
  packages <- as.character(pack[,1])

  for(i in packages){
    library(i, character.only = T)
  }
  rm(i)


  #-----------------

  if("tidyverse" %in% all_packages$Package==FALSE)
    pak::pak("tidyverse/tidyverse")
  library(tidyverse)


  #-----------------
  BiocManager_pack <- data.frame(c("DOSE","clusterProfiler","do","enrichplot",
                                   "pathview","BiocParallel","GO.db","KEGGREST",
                                   "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"))
  # human: "org.Hs.eg.db"
  # mouse: "org.Mn.eg.db"
  # rat: "org.Rn.eg.db"

  BiocManager_pack$type <- BiocManager_pack[,1] %in% all_packages$Package

  for (i in 1:nrow(BiocManager_pack)){
    if(BiocManager_pack[i,2]==FALSE)
      BiocManager::install(BiocManager_pack[i,1],update = F,ask = F)
  }

  # 批量library
  Bio_packages <- as.character(BiocManager_pack[,1])
  for(i in Bio_packages){
    library(i, character.only = T)
  }
  rm(i)


  #--------------------------------------
  #创建KEGG本地数据库
 # if("KEGG.db" %in% all_packages$Package == F)
 # {install_github("YuLab-SMU/createKEGGdb")
 #   library(createKEGGdb)

    #选择创建几个常见物种的kegg注释包: 人"hsa"，小鼠"mmu",大鼠"rno";
  #  kegg_db <-c( "hsa", "mmu", "rno")
  #  createKEGGdb::create_kegg_db(kegg_db)

    #安装这个包(默认的包的路径在当前工作目录，根据实际情况修改路径)
  #  install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")

    #载入自己创建的KEGG.db包；
  #  library(KEGG.db)

  #  file.remove("KEGG.db_1.0.tar.gz")
    #使用本地数据（KEGG.db）进行富集分析
    # 在 enrichKEGG ( use_internal_data= T)
 # }


  #--------------------------------------

  spe <- c("human","mouse","rat")
  species <- tolower(species)

  if(species %in% spe ==FALSE)
  {spec_txt <- paste("The species that you writed was",species,", which was not included in the species of this R package (human, mouse, or rat)!")
  stop(spec_txt)}

  #---------------------------------------

  dir.file <- dplyr::case_when ( tolower(species)== "human" ~ "analysis results (human)",
                                 tolower(species)== "mouse" ~ "analysis results (mouse)",
                                 tolower(species)== "rat" ~ "analysis results (rat)",
                                 TRUE ~ "analysis results"
  )
  if (dir.exists(dir.file)==FALSE)
    dir.create(dir.file)

  group <- paste0("(",group1," VS ",group2,")")
  gene_df_0 <- read.xlsx(gene_file)

  colnames(gene_df_0)[1] <- "gene_id"

  gene_df <- data.frame(distinct(gene_df_0, gene_id, .keep_all = TRUE))

  if(ncol(gene_df)>=2)
  {log <- colnames(gene_df)[2]
  log <- tolower(log)
  log <- stringr::str_sub(log,1,3)
  if(log=="log")
  {colnames(gene_df)[2] <- "log2FC"
  up_gene <- dplyr::filter(gene_df,log2FC>0)
  down_gene <- dplyr::filter(gene_df,log2FC<0)
  up_name <- paste(species,"up_genes_data.xlsx")
  up_name <-paste0(dir.file,"/",up_name)
  write.xlsx(up_gene,up_name)
  down_name <- paste(species,"down_genes_data.xlsx")
  down_name <-paste0(dir.file,"/",down_name)
  write.xlsx(down_gene,down_name)} } else
    print("The gene data did not provide the log2FC.")

  #--------------------------------------------

  gene_n <- dplyr::filter(gene_df,tolower(str_sub(gene_df$gene_id,1,1)) %in% 0:9) # 数字开头的基因
  gene_L <- dplyr::filter(gene_df,tolower(str_sub(gene_df$gene_id,1,1)) %in% letters) # 字母开头的基因

  # hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
  gene_human <- c(gene_n$gene_id,toupper(gene_L$gene_id))

  # mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
  gene_mouse <- c( gene_n$gene_id, str_to_title(tolower(gene_L$gene_id)) )  # str_to_title()首字母大写
  gene_rat <- c( gene_n$gene_id, str_to_title(tolower(gene_L$gene_id)) )

  if(tolower(species)== "human")
  {gene_ENTREZID<-bitr(gene_human,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=org.Hs.eg.db)
  gene_ENTREZID <- na.omit(gene_ENTREZID)
  gene_ENT_h <- paste(species,"gene_ENTREZID (all)_data.xlsx")
  gene_ENT_h <-paste0(dir.file,"/",gene_ENT_h)
  write.xlsx(gene_ENTREZID,gene_ENT_h)
  kegg_gene_df <- enrichKEGG(gene_ENTREZID$ENTREZID, organism = 'hsa',
                             keyType = 'kegg',
                             pvalueCutoff = 0.05,
                             pAdjustMethod = 'BH',
                             qvalueCutoff = 0.2,
                             minGSSize = 3,
                             maxGSSize = 3500,
                             use_internal_data = F) # 是否使用要地kegg数据库，一般是 F

  pathways_geneID <- kegg_gene_df@result

  symbol <- setReadable(kegg_gene_df, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  pathways_geneSymbol <- symbol@result
  kegg_all_pathways <- cbind(pathways_geneID[,1:8],pathways_geneSymbol[,8],pathways_geneID[,9])} else

  {if(tolower(species)== "mouse"){
    gene_ENTREZID<-bitr(gene_mouse,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=org.Mm.eg.db)
    gene_ENTREZID <- na.omit(gene_ENTREZID)
    gene_ENT_m <- paste(species,"gene_ENTREZID (all)_data.xlsx")
    gene_ENT_m <-paste0(dir.file,"/",gene_ENT_m)
    write.xlsx(gene_ENTREZID,gene_ENT_m)

    kegg_gene_df <- enrichKEGG(gene_ENTREZID$ENTREZID, organism ='mmu',
                               #universe,
                               keyType = 'kegg',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',
                               qvalueCutoff = 0.2,
                               minGSSize = 3,
                               maxGSSize = 3500,
                               use_internal_data = F)

    pathways_geneID <- kegg_gene_df@result

    symbol <- setReadable(kegg_gene_df, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    pathways_geneSymbol <- symbol@result
    kegg_all_pathways <- cbind(pathways_geneID[,1:8],pathways_geneSymbol[,8],pathways_geneID[,9])} else

    {if(tolower(species)== "rat"){
      gene_ENTREZID<-bitr(gene_rat,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=org.Rn.eg.db)
      gene_ENTREZID <- na.omit(gene_ENTREZID)
      gene_ENT_r <- paste(species,"gene_ENTREZID (all)_data.xlsx")
      gene_ENT_r <-paste0(dir.file,"/",gene_ENT_r)
      write.xlsx(gene_ENTREZID,gene_ENT_r)
      kegg_gene_df <- enrichKEGG(gene_ENTREZID$ENTREZID, organism = 'rno',
                                 keyType = 'kegg',
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = 'BH',
                                 qvalueCutoff = 0.2,
                                 minGSSize = 3,
                                 maxGSSize = 3500,
                                 use_internal_data = F)
      pathways_geneID <- kegg_gene_df@result
      symbol <- setReadable(kegg_gene_df, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
      pathways_geneSymbol <- symbol@result
      kegg_all_pathways <- cbind(pathways_geneID[,1:8],pathways_geneSymbol[,8],pathways_geneID[,9])} else
        print("The species is error, please check it!")
    }
  }

  #--------------------------------------------

  if(tolower(species)== "human")
  kegg_all_pathways <- tidyr::separate(kegg_all_pathways,Description,"Description","- Homo",remove = T)
  
  if(tolower(species)== "mouse")
    kegg_all_pathways <- tidyr::separate(kegg_all_pathways,Description,"Description","- Mus",remove = T)
  
  if(tolower(species)== "rat")
    kegg_all_pathways <- tidyr::separate(kegg_all_pathways,Description,"Description","- Rat",remove = T)
  
 
  colnames(kegg_all_pathways)[c(1,9,10)] <- c("pathwayID","geneSymbol","Count")

  #把分数字符串变小数
  Gene_Ratio_all <- as.data.frame(apply(str_split(kegg_all_pathways$GeneRatio,"/",simplify = T),2,as.numeric))
  Gene_Ratio_all <- Gene_Ratio_all[,1]/Gene_Ratio_all[,2]

  kegg_all_pathways$GeneRatio <- Gene_Ratio_all

  kegg_all_pathways <- dplyr::arrange(kegg_all_pathways,desc(GeneRatio))

  file_path_name <- paste(species,"KEGG pathways analysis (all)_data.xlsx")
  file_path_name <-paste0(dir.file,"/",file_path_name)

  write.xlsx(kegg_all_pathways,file_path_name)

  row_n <- nrow(kegg_all_pathways)

  if(row_n<30)
    title_gene_text <- paste("The genes enriched pathways","(",group1,"VS",group2,")") else
      title_gene_text <- paste("TOP 30 genes enriched pathways","(",group1,"VS",group2,")")

  title_size <- case_when(nrow(kegg_gene_df)>30 ~12,
                          nrow(kegg_gene_df)>20 ~12,
                          TRUE ~11)

  xy_size <- case_when(nrow(kegg_gene_df)>30 ~7,
                       nrow(kegg_gene_df)>20 ~8,
                       TRUE ~9)

  kegg_mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

  kegg_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=xy_size))+
    theme(legend.text=element_text(face="bold",color="black",size=xy_size))

  if(row_n<30)
    path_n <- row_n else
      path_n <- 30

  kegg_df <- kegg_all_pathways[1:path_n,]


  #--------------kegg all pathways------------------------#

  kegg_pathways <- ggplot(kegg_df)+
    geom_point(aes(x=GeneRatio,
                   y=fct_reorder(Description,GeneRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
    labs(x = 'GeneRatio', y = '',title=title_gene_text)+
    kegg_mytheme+kegg_xytheme

  kegg_pathways

  path_name <- paste(species,"KEGG pathways analysis (all).png")
  path_name <-paste0(dir.file,"/",path_name)

  ggsave(path_name, kegg_pathways,width=1200, height =1000, dpi=150,units = "px")


  #------------kegg metabolism pathways----------------#

  #筛选出含有"metabol"字段的行
  kegg_meta <- dplyr::filter(kegg_all_pathways,grepl("metabol",kegg_all_pathways$Description,ignore.case = T))

  if(nrow(kegg_meta)>0)
  {kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

  title_meta_text <- paste0("KEGG metabolism pathways","(",group1," VS ",group2,")")

  kegg_meta_df <- dplyr::filter(kegg_meta, trimws(tolower(kegg_meta$Description)) !="metabolic pathways")

  kegg_meta_pathways <- ggplot(kegg_meta_df)+
    geom_point(aes(x=GeneRatio,y=reorder(Description,GeneRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title=title_meta_text)+
    kegg_mytheme+kegg_xytheme

  kegg_meta_pathways

  meta_path_name <- paste(species,"metabolism pathways analysis (all).png")
  meta_path_name <-paste0(dir.file,"/",meta_path_name)

  ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

  
  kegg_meta_file <- dplyr::arrange(kegg_meta_df,desc(GeneRatio))

  meta_name <- paste(species,"metabolism pathways analysis (all)_data.xlsx")
  meta_name <-paste0(dir.file,"/",meta_name)

  write.xlsx(kegg_meta_file,meta_name)} else
    print("There is no metabolism pathway!")


  #------------kegg signaling pathways----------------#

  #筛选出含有"signal"字段的行
  kegg_signal <- dplyr::filter(kegg_all_pathways,grepl("signal",kegg_all_pathways$Description,ignore.case = T))
  
  if(nrow(kegg_signal)>0)
  {kegg_signal$Description <- factor(kegg_signal$Description,levels=kegg_signal$Description)

  title_signal_text <- paste0("KEGG signaling pathways"," (",group1," VS ",group2,")")

  kegg_signal_pathways <- ggplot(kegg_signal)+
    geom_point(aes(x=GeneRatio,y=reorder(Description,GeneRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title=title_signal_text)+
    kegg_mytheme+kegg_xytheme

  kegg_signal_pathways

  sign_path_name <- paste(species,"signaling pathways analysis (all).png")
  sign_path_name <-paste0(dir.file,"/",sign_path_name)

  ggsave(sign_path_name, kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")

  
  kegg_signal_file <- dplyr::arrange(kegg_signal,desc(GeneRatio))
  
  sign_name <- paste(species,"signaling pathways analysis (all)_data.xlsx")
  sign_name <-paste0(dir.file,"/",sign_name)

  write.xlsx(kegg_signal_file,sign_name)} else
    print("There is no signaling pathway!")


  
  #=========================================================================#
  
  
  #--------------up-regulated gene pathways------------------------#

  if(ncol(gene_df)>=2){

    if(file.exists(up_name)==TRUE)
    {gene_up_n <- dplyr::filter(up_gene,tolower(str_sub(up_gene$gene_id,1,1)) %in% 0:9) # 数字开头的基因
    gene_up_L <- dplyr::filter(up_gene,tolower(str_sub(up_gene$gene_id,1,1)) %in% letters) # 字母开头的基因

    # hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
    gene_up_human <- c(gene_up_n$gene_id,toupper(gene_up_L$gene_id))

    # mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
    gene_up_mouse <- c( gene_up_n$gene_id, str_to_title(tolower(gene_up_L$gene_id)) )  # str_to_title()首字母大写
    gene_up_rat <- c( gene_up_n$gene_id, str_to_title(tolower(gene_up_L$gene_id)) )

    if(tolower(species)== "human")
    {gene_ENTREZID_up <-bitr(gene_up_human,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
    gene_ENTREZID_up <- na.omit(gene_ENTREZID_up)
    gene_h_u <- paste(species,"gene_ENTREZID (up)_data.xlsx")
    gene_h_u <-paste0(dir.file,"/",gene_h_u)
    write.xlsx(gene_ENTREZID_up,gene_h_u)
    kegg_gene_up <- enrichKEGG(gene_ENTREZID_up$ENTREZID, organism = 'hsa',
                               keyType = 'kegg',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',
                               qvalueCutoff = 0.2,
                               minGSSize = 3,
                               maxGSSize = 3500,
                               use_internal_data = F)
    pathways_geneID_up <- kegg_gene_up@result
    symbol_up <- setReadable(kegg_gene_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    pathways_geneSymbol_up <- symbol_up@result
    kegg_up_pathways <- cbind(pathways_geneID_up[,1:8],pathways_geneSymbol_up[,8],pathways_geneID_up[,9])} else

    {if(tolower(species)== "mouse")
    {gene_ENTREZID_up <-bitr(gene_up_mouse,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
    gene_ENTREZID_up <- na.omit(gene_ENTREZID_up)
    gene_m_u <- paste(species,"gene_ENTREZID (up)_data.xlsx")
    gene_m_u <-paste0(dir.file,"/",gene_m_u)
    write.xlsx(gene_ENTREZID_up,gene_m_u)
    kegg_gene_up <- enrichKEGG(gene_ENTREZID_up$ENTREZID, organism = 'mmu',
                               keyType = 'kegg',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',
                               qvalueCutoff = 0.2,
                               minGSSize = 3,
                               maxGSSize = 3500,
                               use_internal_data = F)
    pathways_geneID_up <- kegg_gene_up@result
    symbol_up <- setReadable(kegg_gene_up, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    pathways_geneSymbol_up <- symbol_up@result
    kegg_up_pathways <- cbind(pathways_geneID_up[,1:8],pathways_geneSymbol_up[,8],pathways_geneID_up[,9])} else

    {if(tolower(species)== "rat")
    {gene_ENTREZID_up <-bitr(gene_up_rat,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Rn.eg.db")
    gene_ENTREZID_up <- na.omit(gene_ENTREZID_up)
    gene_r_u <- paste(species,"gene_ENTREZID (up)_data.xlsx")
    gene_r_u <-paste0(dir.file,"/",gene_r_u)
    write.xlsx(gene_ENTREZID_up,gene_r_u)
    kegg_gene_up <- enrichKEGG(gene_ENTREZID_up$ENTREZID, organism = 'rno',
                               keyType = 'kegg',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',
                               qvalueCutoff = 0.2,
                               minGSSize = 3,
                               maxGSSize = 3500,
                               use_internal_data = F)
    pathways_geneID_up <- kegg_gene_up@result
    symbol_up <- setReadable(kegg_gene_up, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
    pathways_geneSymbol_up <- symbol_up@result
    kegg_up_pathways <- cbind(pathways_geneID_up[,1:8],pathways_geneSymbol_up[,8],pathways_geneID_up[,9])} else
      print("The species is error, please check it!")
    }
    }

    if(tolower(species)== "human")
      kegg_up_pathways <- tidyr::separate(kegg_up_pathways,Description,"Description","- Homo",remove = T)
    
    if(tolower(species)== "mouse")
      kegg_up_pathways <- tidyr::separate(kegg_up_pathways,Description,"Description","- Mus",remove = T)
    
    if(tolower(species)== "rat")
      kegg_up_pathways <- tidyr::separate(kegg_up_pathways,Description,"Description","- Rat",remove = T)
    
    colnames(kegg_up_pathways)[c(1,9,10)] <- c("pathwayID","geneSymbol","Count")

    if(nrow(kegg_up_pathways)>1)
    {Gene_Ratio_up <- data.frame(apply(str_split(kegg_up_pathways$GeneRatio,"/",simplify = T),2,as.numeric))
    Gene_Ratio_up <- Gene_Ratio_up[,1]/Gene_Ratio_up[,2]}else
    {Gene_Ratio_up <-apply(str_split(kegg_up_pathways$GeneRatio,"/",simplify = T),2,as.numeric)
    Gene_Ratio_up <- Gene_Ratio_up[1]/Gene_Ratio_up[2]}

    kegg_up_pathways$GeneRatio <- Gene_Ratio_up

    kegg_up_pathways <- dplyr::arrange(kegg_up_pathways,desc(GeneRatio))

    file_up_name <- paste(species,"up_gene enriched pathways (up)_data.xlsx")
    file_up_name <-paste0(dir.file,"/",file_up_name)

    write.xlsx(kegg_up_pathways,file_up_name)

    row_n_up <- nrow(kegg_up_pathways)

    if(row_n_up<30)
      title_up_text <- paste("The up_genes enriched pathways","(",group1,"VS",group2,")") else
        title_up_text <- paste("TOP 30 up_genes enriched pathways","(",group1,"VS",group2,")")

    title_size_up <- case_when(nrow(kegg_gene_up)>30 ~12,
                               nrow(kegg_gene_up)>20 ~12,
                               TRUE ~11)
    xy_size_up <- case_when(nrow(kegg_gene_up)>30 ~9.5,
                            nrow(kegg_gene_up)>20 ~10.5,
                            TRUE ~11)
    kegg_mytheme_up <-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size_up),
            panel.grid = element_blank(),
            panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
      theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

    kegg_xytheme_up <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size_up,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=xy_size_up))+
      theme(legend.text=element_text(face="bold",color="black",size=xy_size_up))

    if(row_n_up<30)
      path_n_up <- row_n_up else
        path_n_up <- 30

    kegg_up <- kegg_up_pathways[1:path_n_up,]

    kegg_pathways_up <- ggplot(kegg_up)+
      geom_point(aes(x=GeneRatio,
                     y=fct_reorder(Description,GeneRatio),
                     color=-log10(pvalue),size=Count))+
      scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
      labs(x = 'GeneRatio', y = '',title=title_up_text)+
      kegg_mytheme_up+kegg_xytheme_up

    kegg_pathways_up

    path_name_up <- paste(species,"up_genes enriched pathways (up).png")
    path_name_up <-paste0(dir.file,"/",path_name_up)
    ggsave(path_name_up, kegg_pathways_up,width=1200, height =1000, dpi=150,units = "px")}else
      print("up-regulated gene enriched pathways can not be analyaed without up_genes.")



    #--------------down-regulated gene pathways------------------------#

    if(file.exists(down_name)==TRUE)
    {gene_down_n <- dplyr::filter(down_gene,tolower(str_sub(down_gene$gene_id,1,1)) %in% 0:9) # 数字开头的基因
    gene_down_L <- dplyr::filter(down_gene,tolower(str_sub(down_gene$gene_id,1,1)) %in% letters) # 字母开头的基因

    # hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
    gene_down_human <- c(gene_down_n$gene_id,toupper(gene_down_L$gene_id))

    # mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
    gene_down_mouse <- c( gene_down_n$gene_id, str_to_title(tolower(gene_down_L$gene_id)) )  # str_to_title()首字母大写
    gene_down_rat <- c( gene_down_n$gene_id, str_to_title(tolower(gene_down_L$gene_id)) )

    if(tolower(species)== "human")
    {gene_ENTREZID_down <-bitr(gene_down_human,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
    gene_ENTREZID_down <- na.omit(gene_ENTREZID_down)
    gene_h_d <- paste(species,"gene_ENTREZID (down)_data.xlsx")
    gene_h_d <-paste0(dir.file,"/",gene_h_d)
    write.xlsx(gene_ENTREZID_down,gene_h_d)
    kegg_gene_down <- enrichKEGG(gene_ENTREZID_down$ENTREZID, organism = 'hsa',
                                 keyType = 'kegg',
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = 'BH',
                                 qvalueCutoff = 0.2,
                                 minGSSize = 3,
                                 maxGSSize = 3500,
                                 use_internal_data = F)
    pathways_geneID_down <- kegg_gene_down@result
    symbol_down <- setReadable(kegg_gene_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    pathways_geneSymbol_down <- symbol_down@result
    kegg_down_pathways <- cbind(pathways_geneID_down[,1:8],pathways_geneSymbol_down[,8],pathways_geneID_down[,9])} else

    {if(tolower(species)== "mouse")
    {gene_ENTREZID_down <-bitr(gene_down_mouse,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Mm.eg.db")
    gene_ENTREZID_down <- na.omit(gene_ENTREZID_down)
    gene_m_d <- paste(species,"ENTREZID (down)_data.xlsx")
    gene_m_d <-paste0(dir.file,"/",gene_m_d)
    write.xlsx(gene_ENTREZID_down,gene_m_d)
    kegg_gene_down <- enrichKEGG(gene_ENTREZID_down$ENTREZID, organism = 'mmu',
                                 keyType = 'kegg',
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = 'BH',
                                 qvalueCutoff = 0.2,
                                 minGSSize = 3,
                                 maxGSSize = 3500,
                                 use_internal_data = F)
    pathways_geneID_down <- kegg_gene_down@result
    symbol_down <- setReadable(kegg_gene_down, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    pathways_geneSymbol_down <- symbol_down@result
    kegg_down_pathways <- cbind(pathways_geneID_down[,1:8],pathways_geneSymbol_down[,8],pathways_geneID_down[,9])} else

    {if(tolower(species)== "rat")
    {gene_ENTREZID_down <-bitr(gene_down_rat,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Rn.eg.db")
    gene_ENTREZID_down <- na.omit(gene_ENTREZID_down)
    gene_r_d <- paste(species,"gene_ENTREZID (down)_data.xlsx")
    gene_r_d <-paste0(dir.file,"/",gene_r_d)
    write.xlsx(gene_ENTREZID_down,gene_r_d)
    kegg_gene_down <- enrichKEGG(gene_ENTREZID_down$ENTREZID, organism = 'rno',
                                 keyType = 'kegg',
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = 'BH',
                                 qvalueCutoff = 0.2,
                                 minGSSize = 3,
                                 maxGSSize = 3500,
                                 use_internal_data = F)
    pathways_geneID_down <- kegg_gene_down@result
    symbol_down <- setReadable(kegg_gene_down, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
    pathways_geneSymbol_down <- symbol_down@result
    kegg_down_pathways <- cbind(pathways_geneID_down[,1:8],pathways_geneSymbol_down[,8],pathways_geneID_down[,9])} else
      print("The species is error, please check it!")
    }
    }

    if(tolower(species)== "human")
      kegg_down_pathways <- tidyr::separate(kegg_down_pathways,Description,"Description","- Homo",remove = T)
    
    if(tolower(species)== "mouse")
      kegg_down_pathways <- tidyr::separate(kegg_down_pathways,Description,"Description","- Mus",remove = T)
    
    if(tolower(species)== "rat")
      kegg_down_pathways <- tidyr::separate(kegg_down_pathways,Description,"Description","- Rat",remove = T)
    
    colnames(kegg_down_pathways)[c(1,9,10)] <- c("pathwayID","geneSymbol","Count")

    if(nrow(kegg_down_pathways)>1)
    {Gene_Ratio_down <- data.frame(apply(str_split(kegg_down_pathways$GeneRatio,"/",simplify = T),2,as.numeric))
    Gene_Ratio_down <- Gene_Ratio_down[,1]/Gene_Ratio_down[,2]}else
    {Gene_Ratio_down <-apply(str_split(kegg_down_pathways$GeneRatio,"/",simplify = T),2,as.numeric)
    Gene_Ratio_down <- Gene_Ratio_down[1]/Gene_Ratio_down[2]}

    kegg_down_pathways$GeneRatio <- Gene_Ratio_down

    kegg_down_pathways <- dplyr::arrange(kegg_down_pathways,desc(GeneRatio))
    file_down_name <- paste(species,"down_gene enriched pathways (down)_data.xlsx")
    file_down_name <-paste0(dir.file,"/",file_down_name)

    write.xlsx(kegg_down_pathways,file_down_name)

    row_n_down <- nrow(kegg_down_pathways)

    if(row_n_down<30)
      title_down_text <- paste("The down_genes enriched pathways","(",group1,"VS",group2,")") else
        title_down_text <- paste("TOP 30 down_genes enriched pathways","(",group1,"VS",group2,")")

    title_size_down <- case_when(nrow(kegg_gene_down)>30 ~12,
                                 nrow(kegg_gene_down)>20 ~12,
                                 TRUE ~11)
    xy_size_down <- case_when(nrow(kegg_gene_down)>30 ~9.5,
                              nrow(kegg_gene_down)>20 ~10.5,
                              TRUE ~11)
    kegg_mytheme_down <-theme_bw()+
      theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size_down),
            panel.grid = element_blank(),
            panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
            axis.line = element_blank(),
            axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
            axis.ticks.length = unit(1.5,units = "mm"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
      theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

    kegg_xytheme_down <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size_down,angle =0,hjust=1))+
      theme(axis.text.y = element_text(face="bold",color="black",size=xy_size_down))+
      theme(legend.text=element_text(face="bold",color="black",size=xy_size_down))

    if(row_n_down<30)
      path_n_down <- row_n_down else
        path_n_down <- 30

    kegg_down <- kegg_down_pathways[1:path_n_down,]

    kegg_pathways_down <- ggplot(kegg_down)+
      geom_point(aes(x=GeneRatio,
                     y=fct_reorder(Description,GeneRatio),
                     color=-log10(pvalue),size=Count))+
      scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
      labs(x = 'GeneRatio', y = '',title=title_down_text)+
      kegg_mytheme_down+kegg_xytheme_down

    kegg_pathways_down

    path_name_down <- paste(species,"down_genes enriched pathways (down).png")
    path_name_down <-paste0(dir.file,"/",path_name_down)
    ggsave(path_name_down, kegg_pathways_down,width=1200, height =1000, dpi=150,units = "px")}else
      print("down-regulated gene enriched pathways can not be analyaed without down_genes.")

    
  #-----------------up and down-regulation--------------
  
    if(length(kegg_up)>0 & length(kegg_down)>0){
      
      meta_up_n <- grep("Metabolic pathways",kegg_up$Description,ignore.case = T) %>% as.numeric()
      if(length(meta_up_n)>0){
        kegg_up_10 <- kegg_up[-meta_up_n,]}
      
        kegg_up_10 <- dplyr::arrange(kegg_up_10,pvalue) 
        kegg_up_10 <- mutate(kegg_up_10,"log10pvalue"=-log10(pvalue))
      
      meta_down_n <- grep("Metabolic pathways",kegg_down$Description,ignore.case = T) %>% as.numeric()
      if(length(meta_down_n)>0){
        kegg_down_10 <- kegg_down[-meta_down_n,]}
       
        kegg_down_10 <- dplyr::arrange(kegg_down_10,pvalue) 
        kegg_down_10 <- mutate(kegg_down_10,"log10pvalue"=-log10(pvalue))
        kegg_down_10$log10pvalue <- 0-kegg_down_10$log10pvalue
      
      if(nrow(kegg_up_10)>10)
         k_up <- kegg_up_10[1:10,] else
           k_up <- kegg_up_10
        
        if(nrow(kegg_down_10)>10)
          k_down <- kegg_down_10[1:10,] else
            k_down <- kegg_down_10
        
      kegg_10 <- rbind(k_up,k_down)
      kegg_10$Description <- trimws(kegg_10$Description)
      
      Type <- case_when(kegg_10$log10pvalue>=0 ~"up",
                        kegg_10$log10pvalue<0 ~"down")
      
      kegg_10 <- mutate(kegg_10,Type=Type)
      
      kegg_10$log10p_up <- purrr::map(c(1:nrow(kegg_10)),~{if(kegg_10[.x,ncol(kegg_10)]=="up")
                                                               kegg_10[.x,ncol(kegg_10)+1]=round(kegg_10[.x,ncol(kegg_10)-1],2) else
                                                                 kegg_10[.x,ncol(kegg_10)+1]=""  
                                                             } )
      
      kegg_10$log10p_down <- purrr::map(c(1:nrow(kegg_10)),~{if(kegg_10[.x,ncol(kegg_10)-1]=="down")
                                                               kegg_10[.x,ncol(kegg_10)+1]=0-round(kegg_10[.x,ncol(kegg_10)-2],2) else
                                                               kegg_10[.x,ncol(kegg_10)+1]=""  
                                                             } )
     
      ymax <- round(max(kegg_up_10$log10pvalue)+5,0)
      ymin <- round(min(kegg_down_10$log10pvalue)-5,0)
      
      
      
      p0 <- ggplot(kegg_10,aes(x =reorder(Description,log10pvalue),y =log10pvalue,fill=Type))+
        geom_col()+
        theme_classic()+
        ylim(ymin,ymax)+
        coord_flip()+
        scale_fill_manual(values = c("#2f73bb","#ae4531"))+ # 指定颜色
        geom_segment(aes(y=0,yend=0,x=0,xend=nrow(kegg_10)+0.8)) # 加一条坚线
      
      p0
      
      p0+geom_text(aes(label=format(round(log10pvalue, digits = 2), nsmall = 2)), nudge_y= -1, color="black", size =3)
      
      # 调整主题
      my_theme <- theme(
        legend.position = "none", #去除图例
        plot.title=element_text(hjust = 0.5), #标题居中
        axis.line.y = element_blank(), # 隐去y轴线
        axis.title.y = element_blank(), # 隐去y轴"reorder(pathway,value)"
        axis.ticks.y= element_blank(), # 隐去y轴短线
        axis.text.y = element_blank(), # 隐去y轴"pathway1....."
        axis.text = element_text(colour="black",size=10)
        )+
        theme(plot.title = element_text(size=14,hjust = 0.5))+
        theme(axis.title.x = element_text(colour="black",size=14))+
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
      
      # 添加labels
      label_up <- geom_text(data = kegg_10[which(kegg_10$log10pvalue>=0),], #筛选出value>=0，用which函数
                            aes(x=Description,y=0,label=trimws(Description)), # 添加
                            hjust=1, nudge_y =-0.3,  # hjust=1 是右对齐，nudge_y微调距离
                            size=4
                            )
      
      label_down <- geom_text(data = kegg_10[which(kegg_10$log10pvalue<0),], #筛选出value>=0，用which函数
                              aes(x=Description,y=0,label=trimws(Description)), # 添加
                              hjust=0, nudge_y =0.3,  # hjust=0 是左对齐，nudge_y微调距离
                              size=4
                             )
      
      #添加箭头
      arrow_down <- geom_segment(aes(y=-0.2,yend=ymin, x=nrow(kegg_10)+1,xend=nrow(kegg_10)+1),
                                 arrow=arrow(length = unit(0.2,"cm"),type="closed"),
                                 linewidth=0.5
                                 )
      
      
      arrow_up <- geom_segment(aes(y=0.2,yend=ymax, x=nrow(kegg_10)+1,xend=nrow(kegg_10)+1),
                               arrow=arrow(length = unit(0.2,"cm"),type="closed"),
                               linewidth=0.5
                              )
      
      #在柱上显示数值
      col_value_up <- geom_text(aes(label=log10p_up),
                                nudge_y=1.5, # 和柱顶的距离
                                color="black", size =3)
      
      col_value_down <- geom_text(aes(label=log10p_down),
                                nudge_y=-1.5, # 和柱顶的距离
                                color="black", size =3)
      
      
      p01 <- p0+my_theme+
        label_up+label_down+
        arrow_down + 
        arrow_up+ 
        annotate("text",x=nrow(kegg_10)+1.6,y=ymin+1,label="Down",)+
        annotate("text",x=nrow(kegg_10)+1.6,y=ymax-1,label="Up")+
        ylab("-log10(pvalue)")+
        ggtitle("KEGG pathways enrichment\n (up and down-regulation)")+
        scale_x_discrete(expand=expansion(add=c(0.5,2))) # 边距
      
      p01
      
      p02 <- p01+
            col_value_up+
            col_value_down
      
      p02  
      
      path_10_01 <- paste(species,"KEGG enriched pathways (up and down) 01.png")
      path_10_01 <-paste0(dir.file,"/",path_10_01)
      ggsave(path_10_01,p01,width = 1200,height =1000,units = "px",dpi = 150)
      
      path_10_02 <- paste(species,"KEGG enriched pathways (up and down) 02.png")
      path_10_02 <-paste0(dir.file,"/",path_10_02)
      ggsave(path_10_02,p02,width = 1200,height =1000,units = "px",dpi = 150)
      
      }
  }

  #-------------------------------GO enrichment analysis------------------------------#
  if(tolower(species)== "human")
    Org.Db=org.Hs.eg.db else
    {if(tolower(species)== "mouse")
      Org.Db=org.Mm.eg.db else
      {if(tolower(species)== "rat")
        Org.Db=org.Rn.eg.db else
          print("The species if no found and GO can not be analyzed.")
      }
    }

  go_enrich <- enrichGO(gene_ENTREZID$ENTREZID,
                        OrgDb = Org.Db,
                        keyType = 'ENTREZID',
                        ont='ALL',pAdjustMethod = 'BH',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)

  go_result <- data.frame(go_enrich)

  go_result$log10pvalue  <- round(-log10(go_result$pvalue),2)

  go_result <- go_result[order(go_result$ONTOLOGY,-go_result$log10pvalue),]

  GO_file <- paste(species,"GO enrichment (all)_data.xlsx")
  GO_file <-paste0(dir.file,"/",GO_file)
  write.xlsx(go_result,GO_file)

  go_BP_all <- go_result[go_result$ONTOLOGY=="BP",]
  go_BP_p <- go_BP_all[go_BP_all$pvalue<0.05 & go_BP_all$p.adjust<0.05,]

  if(nrow(go_BP_p)>20)
    n_BP=20 else
      n_BP= nrow(go_BP_p)

  go_CC_all <- go_result[go_result$ONTOLOGY=="CC",]
  go_CC_p <- go_CC_all[go_CC_all$pvalue<0.05 & go_CC_all$p.adjust<0.05,]

  if(nrow(go_CC_p)>20)
    n_CC=20 else
      n_CC= nrow(go_CC_p)

  go_MF_all <- go_result[go_result$ONTOLOGY=="MF",]
  go_MF_p <- go_MF_all[go_MF_all$pvalue<0.05 & go_MF_all$p.adjust<0.05,]

  if(nrow(go_MF_p)>20)
    n_MF=20 else
      n_MF= nrow(go_MF_p)

  go_BP <- go_BP_all[c(1:n_BP),]
  go_CC <- go_CC_all[c(1:n_CC),]
  go_MF <- go_MF_all[c(1:n_MF),]

  go_df <- rbind(go_BP,go_CC,go_MF)

  go_df <- go_df[order(go_df$ONTOLOGY,-go_df$log10pvalue),]

  title_size_go <- case_when(nrow(go_df)>50 ~12,
                             TRUE ~11)

  xy_siz_go <- case_when(nrow(go_df)>50 ~9,
                         nrow(go_df)>30 ~10,
                         TRUE ~11)

  go_mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",face="bold",colour ="black",size =title_size_go),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))


  go_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_siz_go,angle =80,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=xy_siz_go))+
    theme(legend.text=element_text(face="bold",color="black",size=xy_siz_go))

  color_text <- c(rep("#20b2aa",nrow(go_BP)),rep("#d2691e",nrow(go_CC)),rep("#6666cc",nrow(go_MF)))
  textcolor_theme <- theme(axis.text.x = element_text(colour = color_text))

  legend_theme01 <- theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.direction = "horizontal",
    legend.position = c(0.5,0.9),
    legend.background = element_blank()
  )

  legend_theme02 <- theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.direction = "horizontal",
    legend.position = "none",
    legend.background = element_blank()
  )


  height_y01 <- max(go_df$log10pvalue)+5

  nn_row <- nrow(go_df) %>% as.numeric()

  go_df$Term <- NA

  for(i in 1:nn_row){
    if(nchar(go_df[i,3])>40)
      go_df[i,12] <- paste(str_sub(go_df[i,3],1,40),"...") else
        go_df[i,12] <- go_df[i,3]
  }

  go_df <- distinct (go_df,Term,.keep_all = T)

  go_df$Term =factor(go_df$Term,levels=go_df$Term)

  legend_text <- c("Biological process(BP)","Cellular component (CC)","Molecular function (MF)")

  GO_plot <- ggplot(go_df, aes(x = Term,y =log10pvalue, fill = ONTOLOGY))

  GO_enrich01 <- GO_plot+ geom_bar(position = "dodge",stat = "identity",width = 0.8)+
    scale_y_continuous(expand = c(0, 0),limits = c(0, height_y01))+
    scale_fill_manual(values = c("#20b2aa", "#d2691e","#6666cc"),label=legend_text)+
    labs(x = '', y = '-log10(pvalue)',title="GO enrichment")+
    go_mytheme+go_xytheme +textcolor_theme+legend_theme01
  # + scale_x_discrete(labels=function(go_enrich) str_wrap(go_enrich, width=50)) #限定x轴字段宽度

  GO_enrich01

  GO_pic <- paste(species,"GO enrichment 01.png")
  GO_pic <-paste0(dir.file,"/",GO_pic)

  ggsave(GO_pic, GO_enrich01,width=1500, height =1000, dpi=150,units = "px")

  GO_f01 <- paste(species,"GO enrichment 01_data.xlsx")
  GO_f01 <-paste0(dir.file,"/",GO_f01)

  write.xlsx(go_df,GO_f01)

  #------------------------GO enrichment 02 --------------------------------------#

  x1=nrow(go_BP)/2
  y1=max(go_df $log10pvalue)+3

  x2=nrow(go_BP)+nrow(go_CC)/2
  y2=max(go_df $log10pvalue)+3

  x3=nrow(go_BP)+nrow(go_CC)+nrow(go_MF)/2
  y3=max(go_df $log10pvalue)+3

  height_y02 <- max(go_df $log10pvalue)+10

  GO_enrich02 <- GO_enrich01+
    theme(panel.border = element_blank(),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_y_continuous(expand = c(0, 0),limits = c(0, height_y02))+

    labs(x = 'GO enrichment',title="",size=18)+

    legend_theme02+

    theme(axis.line.x = element_line(linewidth = 0.8,color = 'black'),
          axis.line.y = element_line(linewidth = 0.8,color = 'white'),
          axis.title.x = element_text(size = 20,color = 'black',face = 'bold'),
          axis.title.y = element_text(size = 14,color = 'black',face = 'bold'),
          plot.title = element_text(size = 20,color = 'black',face = 'bold'),
          legend.text=element_text(size = 14,color = 'black',face = 'bold'))+

    annotate("text", x =x1, y = y1+5, label = "Biological process(BP)",size = 5,
             color="#20b2aa",fontface = "bold")+
    annotate("text", x =x2, y = y2+5, label = "Cellular component (CC)",size = 5,
             color="#d2691e",fontface = "bold")+
    annotate("text", x =x3, y = y3+5, label = "Molecular function (MF)",size = 5,
             color="#6666cc",fontface = "bold")+
    annotation_custom(rectGrob(gp = gpar(col = "#20b2aa",fill = 0, alpha=1,lty = 'dashed',lwd = 1.3)),
                      xmin = 0.5, xmax = 2*x1+0.4, ymin = 0, ymax = y1+1)+
    annotation_custom(rectGrob(gp = gpar(col = "#d2691e",fill = 0, alpha=1,lty = 'dashed',lwd = 1.3)),
                      xmin = 2*x1+0.6, xmax = x2+nrow(go_CC)/2+0.4, ymin = 0, ymax = y2+1)+
    annotation_custom(rectGrob(gp = gpar(col = "#6666cc",fill = 0, alpha=1,lty = 'dashed',lwd = 1.3)),
                      xmin = x3-nrow(go_MF)/2+0.6, xmax = x3++nrow(go_MF)/2+0.5, ymin = 0, ymax = y3+1)


  GO_enrich02

  GO_pic02 <- paste(species,"GO enrichment 02.png")
  GO_pic02 <-paste0(dir.file,"/",GO_pic02)

  ggsave(GO_pic02, GO_enrich02,width=1500, height =1000, dpi=150,units = "px")

  #-------------------------BP analysis---------------------------------#

  go_BP_all <- go_result[go_result$ONTOLOGY=="BP",]
  go_BP_p <- go_BP_all[go_BP_all$pvalue<0.05 & go_BP_all$p.adjust<0.05,]

  if(nrow(go_BP_p)>30)
    n_BP=30 else
      n_BP= nrow(go_BP_p)

  go_BP <- go_BP_all[c(1:n_BP),]

  BP_x_ev <- function(x){eval(parse(text=x))}
  go_BP$Gene_Ratio  <- round(sapply(go_BP$GeneRatio,BP_x_ev),3)

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))


  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=12))+
    theme(legend.text=element_text(face="bold",color="black",size=12))


  go_BP$Term <- NA

  nn_BP <- nrow(go_BP)

  for(i in 1:nn_BP){
    if(nchar(go_BP[i,3])>60)
      go_BP[i,13] <- paste(str_sub(go_BP[i,3],1,40),"...") else
        go_BP[i,13] <- go_BP[i,3]
  }

  go_BP <- dplyr::arrange(go_BP, desc(Gene_Ratio))

  BP_plot <- ggplot(go_BP)+
    geom_point(aes(x=Gene_Ratio,y=fct_reorder(Term,Gene_Ratio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title="GO_BP enrichment")+
    mytheme+xytheme

  BP_plot

  GO_BP <- paste(species,"BP_enrichment.png")
  GO_BP <-paste0(dir.file,"/",GO_BP)
  ggsave(GO_BP, BP_plot,width=1200, height =1000, dpi=150,units = "px")

  GO_BP_ex <- paste(species,"BP_enrichment_data.xlsx")
  GO_BP_ex <-paste0(dir.file,"/",GO_BP_ex)
  write.xlsx(go_BP,GO_BP_ex)

  #-----------------------CC analysis----------------------------------#

  go_CC_all <- go_result[go_result$ONTOLOGY=="CC",]
  go_CC_p <- go_CC_all[go_CC_all$pvalue<0.05 & go_CC_all$p.adjust<0.05,]

  if(nrow(go_CC_p)>30)
    n_CC=30 else
      n_CC= nrow(go_CC_p)

  go_CC <- go_CC_all[c(1:n_CC),]

  CC_x_ev <- function(x){eval(parse(text=x))}
  go_CC$Gene_Ratio  <- round(sapply(go_CC$GeneRatio,CC_x_ev),3)

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))


  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=12))+
    theme(legend.text=element_text(face="bold",color="black",size=12))


  go_CC$Term <- NA

  nn_CC <- nrow(go_CC)

  for(i in 1:nn_CC){
    if(nchar(go_CC[i,3])>60)
      go_CC[i,13] <- paste(str_sub(go_CC[i,3],1,40),"...") else
        go_CC[i,13] <- go_CC[i,3]
  }

  go_CC <- dplyr::arrange(go_CC, desc(Gene_Ratio))

  CC_plot <- ggplot(go_CC)+
    geom_point(aes(x=Gene_Ratio,y=fct_reorder(Term,Gene_Ratio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title="GO_CC enrichment")+
    mytheme+xytheme

  CC_plot

  GO_CC <- paste(species,"CC_enrichment.png")
  GO_CC <-paste0(dir.file,"/",GO_CC)
  ggsave(GO_CC, CC_plot,width=1200, height =1000, dpi=150,units = "px")

  GO_CC_ex <- paste(species,"CC_enrichment_data.xlsx")
  GO_CC_ex <-paste0(dir.file,"/",GO_CC_ex)
  write.xlsx(go_CC,GO_CC_ex)

  #------------------------------MF analysis-------------------------------#

  go_MF_all <- go_result[go_result$ONTOLOGY=="MF",]
  go_MF_p <- go_MF_all[go_MF_all$pvalue<0.05 & go_MF_all$p.adjust<0.05,]

  if(nrow(go_MF_p)>30)
    n_MF=30 else
      n_MF= nrow(go_MF_p)

  go_MF <- go_MF_all[c(1:n_MF),]

  MF_x_ev <- function(x){eval(parse(text=x))}
  go_MF$Gene_Ratio  <- round(sapply(go_MF$GeneRatio,MF_x_ev),3)

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =14),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))


  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=12,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=12))+
    theme(legend.text=element_text(face="bold",color="black",size=12))


  go_MF$Term <- NA

  nn_MF <- nrow(go_MF)

  for(i in 1:nn_MF){
    if(nchar(go_MF[i,3])>60)
      go_MF[i,13] <- paste(str_sub(go_MF[i,3],1,40),"...") else
        go_MF[i,13] <- go_MF[i,3]
  }


  go_MF <- dplyr::arrange(go_MF, desc(Gene_Ratio))

  # Term <- factor(go_MF$Term,ordered = T,levels =go_MF$Term)

  MF_plot <- ggplot(go_MF)+
    geom_point(aes(x=Gene_Ratio,y=fct_reorder(Term,Gene_Ratio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title="GO_MF enrichment")+
    mytheme+xytheme

  MF_plot

  GO_MF <- paste(species,"MF_enrichment.png")
  GO_MF <-paste0(dir.file,"/",GO_MF)
  ggsave(GO_MF, MF_plot,width=1200, height =1000, dpi=150,units = "px")

  GO_MF_ex <- paste(species,"MF_enrichment_data.xlsx")
  GO_MF_ex <-paste0(dir.file,"/",GO_MF_ex)
  write.xlsx(go_MF,GO_MF_ex)

  p <- GO_enrich01+kegg_pathways+plot_layout(ncol=2,nrow=1,widths = c(3,1))

  GO_kegg <- paste(species,"GO enrichment and kegg pathways.png")
  GO_kegg <-paste0(dir.file,"/",GO_kegg)
  ggsave(GO_kegg, p, width=2400, height =1000, dpi=150,units = "px")

  dev.new(width=18,height=7, noRStudioGD = T)

  print("The analysis results can be found in the folder of <analysis result>")

  p

}
