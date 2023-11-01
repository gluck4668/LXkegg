
GO_enrich <- function(kegg_df){

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

gene_ENTREZID <- kegg_df$gene_ENTREZID

dir.file <- kegg_df$dir.file

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

GO_enrich_result <- list(GO_enrich=GO_enrich01)

return(GO_enrich_result)

}

