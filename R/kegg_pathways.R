kegg_pathways <- function(kegg_df){

#--------------kegg all pathways------------------------#

kegg_mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =kegg_df$title_size),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

kegg_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=kegg_df$xy_size,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=kegg_df$xy_size))+
    theme(legend.text=element_text(face="bold",color="black",size=kegg_df$xy_size))

kegg_pathways <- ggplot(kegg_df$kegg_df)+
    geom_point(aes(x=GeneRatio,
                   y=fct_reorder(Description,GeneRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
    labs(x = 'GeneRatio', y = '',title=kegg_df$title_gene_text)+
    kegg_mytheme+kegg_xytheme

kegg_pathways

  path_name <- paste(species,"KEGG pathways analysis (all).png")
  dir.file <- kegg_df$dir.file
  path_name <-paste0(dir.file,"/",path_name)
  ggsave(path_name, kegg_pathways,width=1200, height =1000, dpi=150,units = "px")

#------------kegg metabolism pathways----------------#
  #筛选出含有"metabol"字段的行
  kegg_all_pathways <- kegg_df$kegg_all_pathways

  if(any(grepl("category",colnames(kegg_all_pathways),ignore.case = T))) #判断是否有category
    kegg_meta <- dplyr::filter(kegg_all_pathways,grepl("metabol",kegg_all_pathways$category,ignore.case = T)) else
     kegg_meta <- dplyr::filter(kegg_all_pathways,grepl("metabol",kegg_all_pathways$Description,ignore.case = T))

  if(nrow(kegg_meta)>0)
  {
  kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)
  kegg_meta_df <- dplyr::filter(kegg_meta, trimws(tolower(kegg_meta$Description)) !="metabolic pathways")

  kegg_meta_file <- dplyr::arrange(kegg_meta,desc(GeneRatio))
  meta_name <- paste(species,"metabolism pathways analysis_data.xlsx")
  meta_name <-paste0(dir.file,"/",meta_name)
  write.xlsx(kegg_meta_file,meta_name)

  if(nrow(kegg_meta_df)>=30)
     {n_meta <- 30
     title_meta_text <- paste0("Top 30 KEGG metabolism pathways","(",group1," VS ",group2,")")
     } else
      {n_meta <- nrow(kegg_meta_df)
      title_meta_text <- paste0("KEGG metabolism pathways","(",group1," VS ",group2,")")
      }

kegg_meta_pathways <- ggplot(kegg_meta_df[1:n_meta,])+
    geom_point(aes(x=GeneRatio,y=reorder(Description,GeneRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title=title_meta_text)+
    kegg_mytheme+kegg_xytheme

  kegg_meta_pathways

  meta_path_name <- paste(species,"metabolism pathways analysis.png")
  meta_path_name <-paste0(dir.file,"/",meta_path_name)
  ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

  } else
    print("There is no metabolism pathway!")


  #------------kegg signaling pathways----------------#

  #筛选出含有"signal"字段的行
  kegg_signal <- dplyr::filter(kegg_all_pathways,grepl("signal",kegg_all_pathways$Description,ignore.case = T))

  if(nrow(kegg_signal)>0)
  {
  kegg_signal$Description <- factor(kegg_signal$Description,levels=kegg_signal$Description)

  kegg_signal_file <- dplyr::arrange(kegg_signal,desc(GeneRatio))
  sign_name <- paste(species,"signaling pathways analysis_data.xlsx")
  sign_name <-paste0(dir.file,"/",sign_name)
  write.xlsx(kegg_signal_file,sign_name)

  if(nrow(kegg_signal)>=30)
   {n_sign <- 30
   title_signal_text <- paste0("Top 30 KEGG signaling pathways"," (",group1," VS ",group2,")")
   } else
   {n_sign <- nrow(kegg_signal)
   title_signal_text <- paste0("KEGG signaling pathways"," (",group1," VS ",group2,")")
    }

kegg_signal_pathways <- ggplot(kegg_signal[1:n_sign,])+
    geom_point(aes(x=GeneRatio,y=reorder(Description,GeneRatio),
                   color=-log10(pvalue),size=Count))+
    scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
    labs(x = 'GeneRatio', y = '',title=title_signal_text)+
    kegg_mytheme+kegg_xytheme

  kegg_signal_pathways

  sign_path_name <- paste(species,"signaling pathways analysis.png")
  sign_path_name <-paste0(dir.file,"/",sign_path_name)
  ggsave(sign_path_name, kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")

  } else
    print("There is no signaling pathway!")

#--------------up-regulated gene pathways------------------------#
  if(length(kegg_df$up_gene)>0){

    gene_ENTREZID_up <-bitr(kegg_df$up_gene$gene_id,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=Org_Db)
    gene_ENTREZID_up <- na.omit(gene_ENTREZID_up)
    gene_h_u <- paste(species,"gene_ENTREZID (up)_data.xlsx")
    gene_h_u <-paste0(dir.file,"/",gene_h_u)
    write.xlsx(gene_ENTREZID_up,gene_h_u)

  kegg_gene_up <- enrichKEGG(gene_ENTREZID_up$ENTREZID, organism = organ,
                               keyType = 'kegg',
                               pvalueCutoff = 0.05,
                               pAdjustMethod = 'BH',
                               qvalueCutoff = 0.2,
                               minGSSize = 3,
                               maxGSSize = 3500,
                               use_internal_data = F)
    pathways_geneID_up <- kegg_gene_up@result
    symbol_up <- setReadable(kegg_gene_up, OrgDb = Org_Db, keyType="ENTREZID")
    kegg_up_pathways <- symbol_up@result

    if(tolower(species)== "human")
      kegg_up_pathways <- tidyr::separate(kegg_up_pathways,Description,"Description","- Homo",remove = T)

    if(tolower(species)== "mouse")
      kegg_up_pathways <- tidyr::separate(kegg_up_pathways,Description,"Description","- Mus",remove = T)

    if(tolower(species)== "rat")
      kegg_up_pathways <- tidyr::separate(kegg_up_pathways,Description,"Description","- Rat",remove = T)

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
        title_up_text <- paste("Top 30 up_genes enriched pathways","(",group1,"VS",group2,")")

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
    ggsave(path_name_up, kegg_pathways_up,width=1200, height =1000, dpi=150,units = "px")

    }else

    {kegg_up <- NULL
      print("up-regulated gene enriched pathways can not be analyaed without up_genes.")
     }


#--------------down-regulated gene pathways------------------------#

  if(length(kegg_df$down_gene)>0)
    {gene_ENTREZID_down <-bitr(kegg_df$down_gene$gene_id,fromType="SYMBOL",
                              toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb=Org_Db)
    gene_ENTREZID_down <- na.omit(gene_ENTREZID_down)
    gene_h_d <- paste(species,"gene_ENTREZID (down)_data.xlsx")
    gene_h_d <-paste0(dir.file,"/",gene_h_d)
    write.xlsx(gene_ENTREZID_down,gene_h_d)
    kegg_gene_down <- enrichKEGG(gene_ENTREZID_down$ENTREZID, organism = organ,
                                 keyType = 'kegg',
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = 'BH',
                                 qvalueCutoff = 0.2,
                                 minGSSize = 3,
                                 maxGSSize = 3500,
                                 use_internal_data = F)
    pathways_geneID_down <- kegg_gene_down@result
    symbol_down <- setReadable(kegg_gene_down, OrgDb = Org_Db, keyType="ENTREZID")
    kegg_down_pathways <- symbol_down@result

    if(tolower(species)== "human")
      kegg_down_pathways <- tidyr::separate(kegg_down_pathways,Description,"Description","- Homo",remove = T)

    if(tolower(species)== "mouse")
      kegg_down_pathways <- tidyr::separate(kegg_down_pathways,Description,"Description","- Mus",remove = T)

    if(tolower(species)== "rat")
      kegg_down_pathways <- tidyr::separate(kegg_down_pathways,Description,"Description","- Rat",remove = T)

if(nrow(kegg_down_pathways)>1)
 str_num <- function(i){eval(str2expression(i))} # 转换为小数点
  kegg_down_pathways$GeneRatio <- sapply(kegg_down_pathways$GeneRatio,str_num)

    kegg_down_pathways <- dplyr::arrange(kegg_down_pathways,desc(GeneRatio))
    file_down_name <- paste(species,"down_gene enriched pathways (down)_data.xlsx")
    file_down_name <-paste0(dir.file,"/",file_down_name)
    write.xlsx(kegg_down_pathways,file_down_name)

    row_n_down <- nrow(kegg_down_pathways)

    if(row_n_down<30)
      title_down_text <- paste("The down_genes enriched pathways","(",group1,"VS",group2,")") else
        title_down_text <- paste("Top 30 down_genes enriched pathways","(",group1,"VS",group2,")")

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
      { kegg_down <- NULL
        print("down-regulated gene enriched pathways can not be analyaed without down_genes.")}

#-----------------up and down-regulation--------------

    if(length(kegg_up)>0 & length(kegg_down)>0){

      meta_up_n <- grep("Metabolic pathways",kegg_up$Description,ignore.case = T) %>% as.numeric()
      if(length(meta_up_n)>0)
        kegg_up_20 <- kegg_up[-meta_up_n,] else
          kegg_up_20 <- kegg_up

      kegg_up_20 <- dplyr::arrange(kegg_up_20,pvalue)
      kegg_up_20 <- mutate(kegg_up_20,"log10pvalue"=-log10(pvalue))

      meta_down_n <- grep("Metabolic pathways",kegg_down$Description,ignore.case = T) %>% as.numeric()
      if(length(meta_down_n)>0)
        kegg_down_20 <- kegg_down[-meta_down_n,] else
          kegg_down_20 <- kegg_down

      kegg_down_20 <- dplyr::arrange(kegg_down_20,pvalue)
      kegg_down_20 <- mutate(kegg_down_20,"log10pvalue"=-log10(pvalue))
      kegg_down_20$log10pvalue <- 0-kegg_down_20$log10pvalue

      if(nrow(kegg_up_20)>=20)
        k_up <- kegg_up_20[1:20,] else
          k_up <- kegg_up_20

      if(nrow(kegg_down_20)>=20)
        k_down <- kegg_down_20[1:20,] else
          k_down <- kegg_down_20

      kegg_20 <- rbind(k_up,k_down)
      kegg_20$Description <- trimws(kegg_20$Description)

      Type <- case_when(kegg_20$log10pvalue>=0 ~"up",
                        kegg_20$log10pvalue<0 ~"down")

      kegg_20 <- mutate(kegg_20,Type=Type)

      kegg_20$log10p_up <- purrr::map(c(1:nrow(kegg_20)),~{if(kegg_20[.x,ncol(kegg_20)]=="up")
        kegg_20[.x,ncol(kegg_20)+1]=round(kegg_20[.x,ncol(kegg_20)-1],2) else
          kegg_20[.x,ncol(kegg_20)+1]=""
      } )

      kegg_20$log10p_down <- purrr::map(c(1:nrow(kegg_20)),~{if(kegg_20[.x,ncol(kegg_20)-1]=="down")
        kegg_20[.x,ncol(kegg_20)+1]=0-round(kegg_20[.x,ncol(kegg_20)-2],2) else
          kegg_20[.x,ncol(kegg_20)+1]=""
      } )

      ymax <- round(max(kegg_up_20$log10pvalue)+5,0)
      ymin <- round(min(kegg_down_20$log10pvalue)-5,0)



  p0 <- ggplot(kegg_20,aes(x =reorder(Description,log10pvalue),y =log10pvalue,fill=Type))+
        geom_col()+
        theme_classic()+
        ylim(ymin,ymax)+
        coord_flip()+
        scale_fill_manual(values = c("#2f73bb","#ae4531"))+ # 指定颜色
        geom_segment(aes(y=0,yend=0,x=0,xend=nrow(kegg_20)+0.8)) # 加一条坚线

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
      label_up <- geom_text(data = kegg_20[which(kegg_20$log10pvalue>=0),], #筛选出value>=0，用which函数
                            aes(x=Description,y=0,label=trimws(Description)), # 添加
                            hjust=1, nudge_y =-0.3,  # hjust=1 是右对齐，nudge_y微调距离
                            size=3 # 字体大小
                            )

      label_down <- geom_text(data = kegg_20[which(kegg_20$log10pvalue<0),], #筛选出value>=0，用which函数
                              aes(x=Description,y=0,label=trimws(Description)), # 添加
                              hjust=0, nudge_y =0.3,  # hjust=0 是左对齐，nudge_y微调距离
                              size=3 # 字体大小
                              )

      #添加箭头
      arrow_down <- geom_segment(aes(y=-0.2,yend=ymin, x=nrow(kegg_20)+1,xend=nrow(kegg_20)+1),
                                 arrow=arrow(length = unit(0.2,"cm"),type="closed"),
                                 linewidth=0.5
                              )


      arrow_up <- geom_segment(aes(y=0.2,yend=ymax, x=nrow(kegg_20)+1,xend=nrow(kegg_20)+1),
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
        annotate("text",x=nrow(kegg_20)+1.6,y=ymin+3,label="Down",)+
        annotate("text",x=nrow(kegg_20)+1.6,y=ymax-3,label="Up")+
        ylab("-log10(pvalue)")+
        ggtitle("KEGG pathways enrichment\n (up and down-regulation)")+
        scale_x_discrete(expand=expansion(add=c(0.5,2))) # 边距

      p01

      p02 <- p01+
        col_value_up+
        col_value_down

      p02

      path_20_01 <- paste(species,"KEGG enriched pathways (up and down) 01.png")
      path_20_01 <-paste0(dir.file,"/",path_20_01)
      ggsave(path_20_01,p01,width = 1200,height =1000,units = "px",dpi = 150)

      path_20_02 <- paste(species,"KEGG enriched pathways (up and down) 02.png")
      path_20_02 <-paste0(dir.file,"/",path_20_02)
      ggsave(path_20_02,p02,width = 1200,height =1000,units = "px",dpi = 150)



    }

  pathways_result <- list(kegg_pathways=kegg_pathways)
  return(pathways_result)

   }



