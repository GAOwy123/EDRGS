rm(list=ls())
library(tidyverse)
setwd(".")
dt<-read_tsv("cellline_NTP.txt")
dt<-dt[,c("cellline","methods","prediction")]



cellline_count_summary<-dt %>% group_by(cellline) %>% summarise(Total_Num=n())
cellline_summary<-dt %>% group_by(cellline,prediction) %>% summarise(count=n())
dt2<-inner_join(cellline_count_summary,cellline_summary)
dt3<-dt2 %>% filter(count>Total_Num/2)

dt4<-data.frame(cellline=dt3$cellline,methods=rep("Final",nrow(dt3)),prediction=dt3$prediction)
dt4<-dt4 %>% dplyr::arrange(prediction)

dt<-dt %>% filter(cellline %in% dt4$cellline)
dt<-rbind(dt4,dt)
cellline_order<-c("TE1", "TE4","TE8", "TE14","TE15",  "HCE4" , "KYSE150",  "KYSE220" , "KYSE270" , "OE19",   
                   "COLO680N", "ECGI10","KYSE140" , "KYSE180","KYSE410",  "KYSE450",  "KYSE50","KYSE510" ,
                  "KYSE520",  "KYSE70","OE21", "TDOTT", "TE11","TE12","TE5","TE6","TE9")
dt$cellline<-factor(dt$cellline,levels=cellline_order)
dt$methods<-factor(dt$methods,levels = c("GDSC1","GDSC2","PRISM","CGP2016","CTRP2","Final"))
dt$prediction<-factor(dt$prediction,levels = c("Low","High"))
dt<-as.data.frame(dt)

mat<-dt %>% spread(cellline,prediction)
df<-reshape2::melt(mat)
pdf("output/14.drug_sensitivity_GDSC/ESCC_cellline/cellline_NTP.pdf",width=8,height = 3)
ggplot() + 
  geom_tile(data=dt, aes(x=cellline, y=methods, fill=factor(prediction)),
            width=1,height=1) + 
  scale_fill_manual(values=c("#AED9E0","#FFA69E"))+#"#46b566","#f08d85"
  theme_bw()+
  theme(panel.background = element_rect(fill="#EAF3DD",colour = "black",size=1.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,color = "black"),
        axis.text.y=element_text(color="black"),
        axis.ticks = element_blank())
dev.off()
