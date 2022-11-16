

#----------------KEGG analysis--------------------------------------#

LXkegg <- function(gene_file,group1,group2){


# list all the packages that have been installed
  all_packages <- data.frame(installed.packages())

# To judge whether a package was installed. If not, it will be installed.
  pack <- data.frame(c("devtools","BiocManager","ggnewscale","tidyverse", "R.utils",
                       "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych","ggplot2",
                       "ggrepel","RColorBrewer", "ggthemes","rticles","grid","patchwork") )

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
                                 "pathview","BiocParallel","org.Hs.eg.db","GO.db"))

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



#---------------------------------------

  group <- paste("(",group1,"VS",group2,")")
  gene_df_0 <- read.xlsx(gene_file)

  colnames(gene_df_0)[1] <- "gene_id"

  gene_df <- data.frame(distinct(gene_df_0, gene_id, .keep_all = TRUE))

  gene_id <- toupper(gene_df[,1])
  #g_frame <- data.frame(gene_id)

  gene_ENTREZID<-bitr(gene_id,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
  gene_ENTREZID <- na.omit(gene_ENTREZID)

  if (dir.exists("analysis results")==FALSE)
      dir.create("analysis results")

  write.xlsx(gene_ENTREZID,"analysis results/gene_ENTREZID.xlsx")

  kegg_gene_df <- enrichKEGG(gene_ENTREZID$ENTREZID, organism = 'hsa',
                           keyType = 'kegg',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = 'BH',
                           qvalueCutoff = 0.2,
                           minGSSize = 3,
                           maxGSSize = 500,
                           use_internal_data = FALSE)

  pathways_geneID <- kegg_gene_df@result

  symbol <- setReadable(kegg_gene_df, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  pathways_geneSymbol <- symbol@result

  kegg_all_pathways <- cbind(pathways_geneID[,1:8],pathways_geneSymbol[,8],pathways_geneID[,9])
  colnames(kegg_all_pathways)[c(1,9,10)] <- c("pathwayID","geneSymbol","Count")
  
  #把分数字符串变小数
  Gene_Ratio <- as.data.frame(apply(str_split(kegg_all_pathways$GeneRatio,"/",simplify = T),2,as.numeric))
  Gene_Ratio <- Gene_Ratio[,1]/Gene_Ratio[,2]
  
  kegg_all_pathways$GeneRatio <- Gene_Ratio
  
  kegg_all_pathways <- dplyr::arrange(kegg_all_pathways,desc(GeneRatio))
  
  write.xlsx(kegg_all_pathways,"analysis results/KEGG pathways analysis.xlsx")

  row_n <- nrow(kegg_all_pathways)

 if(row_n<30)
    title_gene_text <- paste("The genes enriched pathways","(",group1,"VS",group2,")") else
      title_gene_text <- paste("TOP 30 genes enriched pathways","(",group1,"VS",group2,")")

 title_size <- case_when(nrow(kegg_gene_df)>30 ~12,
                          nrow(kegg_gene_df)>20 ~12,
                          TRUE ~11)

 xy_size <- case_when(nrow(kegg_gene_df)>30 ~9,
                       nrow(kegg_gene_df)>20 ~10,
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
  geom_point(aes(x=GeneRatio,y=fct_reorder(Description,Count),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'GeneRatio', y = '',title=title_gene_text)+
  kegg_mytheme+kegg_xytheme

kegg_pathways

ggsave("analysis results/KEGG pathways analysis.png", kegg_pathways,width=1200, height =1000, dpi=150,units = "px")


#------------kegg metabolism pathways----------------#

kegg_meta <- kegg_gene_df

#把Description转换成大写
upper <- toupper(kegg_meta@result$Description)

#筛选出含有"METABOL"字段的行
kegg_meta <- filter(kegg_meta@result,grepl("METABOL",upper))

kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

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

ggsave("analysis results/KEGG metabolism pathways.png", kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

#重新命名
colnames(kegg_all_pathways)[1] <- "Pathways_ID"
colnames(kegg_meta)[1] <- "Pathways_ID"

kegg_meta01 <- kegg_meta[,c(1,10)]
kegg_meta02 <- inner_join(kegg_all_pathways,kegg_meta01,by="Pathways_ID")

#按Gene_Ratio降序排序
kegg_meta02 <- kegg_meta02[order(-kegg_meta02$Gene_Ratio),]
rownames(kegg_meta02) <- c(1:nrow(kegg_meta02))

write.xlsx(kegg_meta02,"analysis results/KEGG metabolism pathways.xlsx")



#------------kegg signaling pathways----------------#

kegg_signal <- kegg_gene_df

#把Description转换成大写
upper <- toupper(kegg_signal@result$Description)

#筛选出含有"METABOL"字段的行
kegg_signal <- filter(kegg_signal@result,grepl("SIGNAL",upper))

kegg_signal$Description <- factor(kegg_signal$Description,levels=kegg_signal$Description)

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

ggsave("analysis results/KEGG signaling pathways.png", kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")

kegg_signal01 <- kegg_signal[order(-kegg_signal$Gene_Ratio),]


#重新命名
#colnames(kegg_all_pathways)[1] <- "Pathways_ID"
colnames(kegg_signal)[1] <- "Pathways_ID"

kegg_signal01 <- kegg_signal[,c(1,10)]
kegg_signal02 <- inner_join(kegg_all_pathways,kegg_signal01,by="Pathways_ID")

#按Gene_Ratio降序排序
kegg_signal02 <- kegg_signal02[order(-kegg_signal02$Gene_Ratio),]
rownames(kegg_signal02) <- c(1:nrow(kegg_signal02))

write.xlsx(kegg_signal02,"analysis results/KEGG signaling pathways.xlsx")



#-------------------------------GO enrichment analysis------------------------------#

  go_enrich <- enrichGO(gene_ENTREZID$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'ENTREZID',
                      ont='ALL',pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

  go_result <- data.frame(go_enrich)

  go_result$log10pvalue  <- round(-log10(go_result$pvalue),2)
  go_result <- go_result[order(go_result$ONTOLOGY,-go_result$log10pvalue),]

  write.xlsx(go_result,"analysis results/GO enrichment (all).xlsx")

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

ggsave("analysis results/GO enrichment 01.png", GO_enrich01,width=1500, height =1000, dpi=150,units = "px")

write.xlsx(go_df,"analysis results/GO enrichment 01.xlsx")

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

              theme(axis.line.x = element_line(size = 0.8,color = 'black'),
                    axis.line.y = element_line(size = 0.8,color = 'white'),
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

  ggsave("analysis results/GO enrichment 02.png", GO_enrich02,width=1500, height =1000, dpi=150,units = "px")

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
    theme(panel.grid =element_line(colour="#dcdcdc",size=0.2,linetype = "dashed"))


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

ggsave("analysis results/BP_enrichment.png", BP_plot,width=1200, height =1000, dpi=150,units = "px")

write.xlsx(go_BP,"analysis results/BP_enrichment.xlsx")

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
    theme(panel.grid =element_line(colour="#dcdcdc",size=0.2,linetype = "dashed"))


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

  ggsave("analysis results/CC_enrichment.png", CC_plot,width=1200, height =1000, dpi=150,units = "px")

  write.xlsx(go_CC,"analysis results/CC_enrichment.xlsx")

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
    theme(panel.grid =element_line(colour="#dcdcdc",size=0.2,linetype = "dashed"))


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

  ggsave("analysis results/MF_enrichment.png", MF_plot,width=1200, height =1000, dpi=150,units = "px")

  write.xlsx(go_MF,"analysis results/MF_enrichment.xlsx")


  p <- GO_enrich01+kegg_pathways+plot_layout(ncol=2,nrow=1,widths = c(3,1))

  ggsave("analysis results/GO enrichment and kegg pathways.png", p, width=2400, height =1000, dpi=150,units = "px")

  dev.new(width=18,height=7, noRStudioGD = T)

  print("The analysis results can be found in the folder of <analysis result>")

  p

}

