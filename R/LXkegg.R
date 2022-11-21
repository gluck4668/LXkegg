

#----------------KEGG analysis--------------------------------------#

LXkegg <- function(gene_file,group1,group2,species){

 #list all the packages that have been installed
  all_packages <- data.frame(installed.packages())

 #To judge whether a package was installed. If not, it will be installed.
  pack <- data.frame(c("devtools","BiocManager","ggnewscale","tidyverse", "R.utils",
                       "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych","ggplot2",
                       "ggrepel","RColorBrewer", "ggthemes","rticles","grid","patchwork","Hmisc") )

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

  

  pack2 <- data.frame(c("tidyverse") )

  # To judge whether a package was included in the all_packages: %in%
  pack$type <- pack2[,1] %in% all_packages$Package

  for (i in 1:nrow(pack2)){
    if(pack[i,2]==FALSE)
      install.packages(pack2[i,1],update = F,ask = F)
  }
  rm(i)

  # 批量library
  packages <- as.character(pack2[,1])

  for(i in packages){
    library(i, character.only = T)
  }
  rm(i)


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
  if("KEGG.db" %in% all_packages$Package == F)
  {install_github("YuLab-SMU/createKEGGdb")
  library(createKEGGdb)
  
  #选择创建几个常见物种的kegg注释包: 人"hsa"，小鼠"mmu",大鼠"rno";
  kegg_db <-c( "hsa", "mmu", "rno")
  createKEGGdb::create_kegg_db(kegg_db)
  
  #安装这个包(默认的包的路径在当前工作目录，根据实际情况修改路径)
  install.packages("KEGG.db_1.0.tar.gz",repos=NULL,type="source")
  
  #载入自己创建的KEGG.db包；
  library(KEGG.db)
  
  file.remove("KEGG.db_1.0.tar.gz")
  #使用本地数据（KEGG.db）进行富集分析
  # 在 enrichKEGG ( use_internal_data= T)
  }
  
  
  #--------------------------------------    

  spe <- c("human","mouse","rat")
  species <- tolower(species)

  if(species %in% spe ==TRUE){

#---------------------------------------

  dir.file <- dplyr::case_when ( tolower(species)== "human" ~ "analysis results (human)",
                                tolower(species)== "mouse" ~ "analysis results (mouse)",
                                tolower(species)== "rat" ~ "analysis results (rat)",
                                TRUE ~ "analysis results"
                               )
 if (dir.exists(dir.file)==FALSE)
    dir.create(dir.file)

  group <- paste("(",group1,"VS",group2,")")
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
      print("The gene data did not provide the log2FoldChange.")

  #--------------------------------------------

  gene_n <- filter(gene_df,tolower(str_sub(gene_df$gene_id,1,1)) %in% 0:9) # 数字开头的基因
  gene_L <- filter(gene_df,tolower(str_sub(gene_df$gene_id,1,1)) %in% letters) # 字母开头的基因
  
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
                           use_internal_data = T)
  
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
                               use_internal_data = T)
      
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
                                use_internal_data = T)
     pathways_geneID <- kegg_gene_df@result
     symbol <- setReadable(kegg_gene_df, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
     pathways_geneSymbol <- symbol@result
     kegg_all_pathways <- cbind(pathways_geneID[,1:8],pathways_geneSymbol[,8],pathways_geneID[,9])} else
     print("The species is error, please check it!")
     }
    }

 #--------------------------------------------

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

 xy_size <- case_when(nrow(kegg_gene_df)>30 ~9.5,
                       nrow(kegg_gene_df)>20 ~10.5,
                       TRUE ~11)

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

kegg_meta <- kegg_gene_df

#把Description转换成大写
upper <- toupper(kegg_meta@result$Description)

#筛选出含有"METABOL"字段的行
kegg_meta <- filter(kegg_meta@result,grepl("METABOL",upper))

if(nrow(kegg_meta)>0)
{kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

#把分数字符串变小数
Gene_Ratio <- as.data.frame(apply(str_split(kegg_meta$GeneRatio,"/",simplify = T),2,as.numeric))
Gene_Ratio <- Gene_Ratio[,1]/Gene_Ratio[,2]

kegg_meta$Gene_Ratio <- Gene_Ratio

title_meta_text <- paste("KEGG metabolism pathways","(",group1,"VS",group2,")")

kegg_meta_pathways <- ggplot(kegg_meta)+
  geom_point(aes(x=Gene_Ratio,y=reorder(Description,Gene_Ratio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
  labs(x = 'GeneRatio', y = '',title=title_meta_text)+
  kegg_mytheme+kegg_xytheme

kegg_meta_pathways

meta_path_name <- paste(species,"metabolism pathways analysis (all).png")
meta_path_name <-paste0(dir.file,"/",meta_path_name)

ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

#重新命名
colnames(kegg_all_pathways)[1] <- "Pathways_ID"
colnames(kegg_meta)[1] <- "Pathways_ID"

kegg_meta01 <- kegg_meta[,c(1,10)]
kegg_meta02 <- inner_join(kegg_all_pathways,kegg_meta01,by="Pathways_ID")

#按Gene_Ratio降序排序
kegg_meta02 <- kegg_meta02[order(-kegg_meta02$Gene_Ratio),]
rownames(kegg_meta02) <- c(1:nrow(kegg_meta02))

kegg_meta02 <- dplyr::arrange(kegg_meta02,desc(Gene_Ratio))

meta_name <- paste(species,"metabolism pathways analysis (all)_data.xlsx")
meta_name <-paste0(dir.file,"/",meta_name)

write.xlsx(kegg_meta02,meta_name)} else
  print("There is no metabolism pathway!")


#------------kegg signaling pathways----------------#

kegg_signal <- kegg_gene_df

#把Description转换成大写
upper <- toupper(kegg_signal@result$Description)

#筛选出含有"METABOL"字段的行
kegg_signal <- filter(kegg_signal@result,grepl("SIGNAL",upper))

if(nrow(kegg_signal)>0)
{kegg_signal$Description <- factor(kegg_signal$Description,levels=kegg_signal$Description)

#把分数字符串变小数
Gene_Ratio <- as.data.frame(apply(str_split(kegg_signal$GeneRatio,"/",simplify = T),2,as.numeric))
Gene_Ratio <- Gene_Ratio[,1]/Gene_Ratio[,2]

kegg_signal$Gene_Ratio <- Gene_Ratio

title_signal_text <- paste("KEGG signaling pathways","(",group1,"VS",group2,")")

kegg_signal_pathways <- ggplot(kegg_signal)+
  geom_point(aes(x=Gene_Ratio,y=reorder(Description,Gene_Ratio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
  labs(x = 'GeneRatio', y = '',title=title_signal_text)+
  kegg_mytheme+kegg_xytheme

kegg_signal_pathways

sign_path_name <- paste(species,"signaling pathways analysis (all).png")
sign_path_name <-paste0(dir.file,"/",sign_path_name)

ggsave(sign_path_name, kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")

kegg_signal01 <- kegg_signal[order(-kegg_signal$Gene_Ratio),]


#重新命名
#colnames(kegg_all_pathways)[1] <- "Pathways_ID"
colnames(kegg_signal)[1] <- "Pathways_ID"

kegg_signal01 <- kegg_signal[,c(1,10)]
kegg_signal02 <- inner_join(kegg_all_pathways,kegg_signal01,by="Pathways_ID")

#按Gene_Ratio降序排序
kegg_signal02 <- kegg_signal02[order(-kegg_signal02$Gene_Ratio),]
rownames(kegg_signal02) <- c(1:nrow(kegg_signal02))

kegg_signal02 <- dplyr::arrange(kegg_signal02,desc(Gene_Ratio))

sign_name <- paste(species,"signaling pathways analysis (all)_data.xlsx")
sign_name <-paste0(dir.file,"/",sign_name)

write.xlsx(kegg_signal02,sign_name)} else
  print("There is no signaling pathway!")


#--------------up-regulated gene pathways------------------------#

if(ncol(gene_df)>=2){

if(file.exists(up_name)==TRUE)
  {gene_up_n <- filter(up_gene,tolower(str_sub(up_gene$gene_id,1,1)) %in% 0:9) # 数字开头的基因
   gene_up_L <- filter(up_gene,tolower(str_sub(up_gene$gene_id,1,1)) %in% letters) # 字母开头的基因
  
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
                               use_internal_data = T)
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
                                 use_internal_data = T)
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
                                   use_internal_data = T)
        pathways_geneID_up <- kegg_gene_up@result
        symbol_up <- setReadable(kegg_gene_up, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
        pathways_geneSymbol_up <- symbol_up@result
        kegg_up_pathways <- cbind(pathways_geneID_up[,1:8],pathways_geneSymbol_up[,8],pathways_geneID_up[,9])} else
          print("The species is error, please check it!")
         }
    }

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
{gene_down_n <- filter(down_gene,tolower(str_sub(down_gene$gene_id,1,1)) %in% 0:9) # 数字开头的基因
gene_down_L <- filter(down_gene,tolower(str_sub(down_gene$gene_id,1,1)) %in% letters) # 字母开头的基因

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
                           use_internal_data = T)
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
                           use_internal_data = T)
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
                           use_internal_data = T)
pathways_geneID_down <- kegg_gene_down@result
symbol_down <- setReadable(kegg_gene_down, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
pathways_geneSymbol_down <- symbol_down@result
kegg_down_pathways <- cbind(pathways_geneID_down[,1:8],pathways_geneSymbol_down[,8],pathways_geneID_down[,9])} else
  print("The species is error, please check it!")
}
}

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

  p} else
    print("The species is error, please check it!")
}


