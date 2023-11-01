
kegg_data_processing <- function(gene_file,group1,group2,species) {
  
  #list all the packages that have been installed
  all_packages <- data.frame(installed.packages())
  
  #To judge whether a package was installed. If not, it will be installed.
  pack <- c("devtools","BiocManager","ggnewscale","R.utils","roxygen2","xfun", 
            "ggsci","openxlsx","dplyr","psych","ggplot2","ggrepel","purrr",
            "RColorBrewer", "ggthemes","rticles","grid","patchwork","Hmisc","pak")
  
  # To judge whether a package was included in the all_packages: %in%
  not_install_pack <- pack[!pack %in% all_packages$Package]
  
  # To install the packages using the sapply function
  if (length(not_install_pack)>0){
    install_fun <- function(x){install.packages(x)}
    sapply(not_install_pack,install_fun)
    rm(install_fun)
  }
  
  #-----------------
  if(!"tidyverse" %in% all_packages$Package)
    pak::pak("tidyverse/tidyverse")
  library(tidyverse)
  
  #-----------------
  Bio_pack <-c("DOSE","clusterProfiler","do","enrichplot",
               "pathview","BiocParallel","GO.db","KEGGREST",
               "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db")
  # human: "org.Hs.eg.db"
  # mouse: "org.Mn.eg.db"
  # rat: "org.Rn.eg.db"
  not_install_bio <- Bio_pack[!Bio_pack %in% all_packages$Package]
  
  if(length(not_install_bio)>0){
    bio_fun <- function(x){BiocManager::install(not_install_bio,update = F,ask = F)}
    sapply(not_install_bio,bio_fun)
    rm(bio_fun)
  }
  
  
  # library packages
  lib_fun <- function(x){library(x,character.only = T)}
  sapply(c(pack,Bio_pack),lib_fun)
  
  rm(all_packages,pack,not_install_pack,
     Bio_pack,not_install_bio,lib_fun)
  
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
  
  if(!species %in% spe)
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
 
  # hsa: all letters of the gene symbol must be UPPER（人类：全部大写）
  if(tolower(species)== "human")
    gene_df$gene_id <- toupper(gene_df$gene_id) else
    
      # mouse/rat: the first letter is UPPER, and the rest of letters must be LOWER（鼠：首大写，其余小写）
      {
       gene_n <- dplyr::filter(gene_df,!is.na(str_extract(gene_df$gene_id,"^\\d+"))) # 数字开头的基因
       gene_n$gene_id <- toupper(gene_n$gene_id)  
       
       id_str <- str_extract(gene_n$gene_id,"[A-Z]+$")
       
       gene_n0 <- gene_n[grep(TRUE,is.na(id_str)),] %>% data.frame()
       colnames(gene_n0)[1] <- "gene_id"
       gene_n1 <- gene_n[grep(FALSE,is.na(id_str)),] %>% data.frame()
       colnames(gene_n1)[1] <- "gene_id"
       
       id_str_1 <- str_extract(gene_n1$gene_id,"[A-Z]+$")
       
       gene_n1$gene_id <- str_replace(gene_n1$gene_id,id_str_1, str_to_title(id_str_1))
  
  
       gene_L <- dplyr::filter(gene_df,!is.na(str_extract(gene_df$gene_id,"^\\D+"))) # 字母开头的基因
       gene_L$gene_id <- str_to_title(gene_L$gene_id) # str_to_title()首字母大写
      
       gene_df <- rbind(gene_n1,gene_n0,gene_L)
       }
   
  
  # 判断是否有log2FC，并筛选出up_gene和down_gene
  if(ncol(gene_df)>=2)
  {log <- grep("log",colnames(gene_df),ignore.case = T)
      if(length(log)>0)
       {colnames(gene_df)[log] <- "log2FC"
        up_gene <- dplyr::filter(gene_df,log2FC>0)
        down_gene <- dplyr::filter(gene_df,log2FC<0)
        up_name <- paste(species,"up_genes_data.xlsx")
        up_name <-paste0(dir.file,"/",up_name)
        write.xlsx(up_gene,up_name)
        down_name <- paste(species,"down_genes_data.xlsx")
        down_name <-paste0(dir.file,"/",down_name)
        write.xlsx(down_gene,down_name)} } else
  {print("The gene data did not provide the log2FC.")
    up_gene <- NULL
    down_gene <- NULL
  }
  
  
  if(tolower(species)== "human")
  {gene_ENTREZID<-bitr(gene_df$gene_id,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=org.Hs.eg.db)
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
    gene_ENTREZID<-bitr(gene_df$gene_id,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=org.Mm.eg.db)
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
      gene_ENTREZID<-bitr(gene_df$gene_id,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=org.Rn.eg.db)
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
  ratio_fun <- function(x){eval(parse(text=x))} # eval(parse(text=x)) 把字符串当做命令运行
  kegg_all_pathways$GeneRatio <- purrr::map(kegg_all_pathways$GeneRatio,ratio_fun) %>% as.numeric()
  
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
  
  
  
  if(row_n<30)
    path_n <- row_n else
      path_n <- 30
  
  kegg_df <- kegg_all_pathways[1:path_n,]
  
  processing_result <- list(title_size=title_size,
                            xy_size=xy_size,
                            title_gene_text=title_gene_text,
                            dir.file=dir.file,
                            up_gene=up_gene,
                            down_gene=down_gene,
                            kegg_df=kegg_df,
                            kegg_all_pathways=kegg_all_pathways,
                            gene_ENTREZID=gene_ENTREZID
                            )
  
  
  
  
  return(processing_result)
  
}
