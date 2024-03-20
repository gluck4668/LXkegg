
df <- read.xlsx("4 GO_data.xlsx")

bb_num <- filter(df,ONTOLOGY=="BP") %>% nrow()
cc_num <- filter(df,ONTOLOGY=="CC") %>% nrow()
mf_num<- filter(df,ONTOLOGY=="MF") %>% nrow()


if(bb_num>=10)
    bb <- 10 else
      bb <- bb_num
if(cc_num>=10)
      cc <- 10 else
        cc <- cc_num
if(mf_num>=10)
       mf <- 10 else
         mf <- mf_num


