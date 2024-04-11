rm(list=ls())
setwd(".")

GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds')) #17419   805
GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds")) # 805 198
GDSC1_Res <- exp(GDSC1_Res) 
library(tidyverse)

ESCC_cellline<-read_tsv("ESCC_cellline.txt",col_names = F)
##+==========================GDSC1 and GDSC2 extract ESCC cell line expression=====================

GDSC_cellline<-read_csv("Cosmic_Cell_listFri Dec  9 02_05_19 2022.csv")
GDSC_cellline$`Cell line Name`<-gsub("-","",GDSC_cellline$`Cell line Name`)
GDSC_cellline<-GDSC_cellline %>% filter(`Cell line Name` %in% ESCC_cellline$X1)
GDSC1_cellline<-GDSC_cellline %>% filter(Datasets=="GDSC1")

ESCC_GDSC1_Expr<-GDSC1_Expr[,paste0("COSMIC_",GDSC1_cellline$`COSMIC ID`)]
GDSC1_Res<-GDSC1_Res[colnames(ESCC_GDSC1_Expr),]
colnames(ESCC_GDSC1_Expr)<-GDSC1_cellline$`Cell line Name`
rownames(GDSC1_Res)<-GDSC1_cellline$`Cell line Name`
###+====================calculate EDDRGS=================================================
###+=============load phenotype

model<-read_tsv("Model.txt")#gene_final.txt
colnames(model)<-c("gene","coef")

exp<-ESCC_GDSC1_Expr[model$gene,]


expr2<-t(exp)
expr_test2<-expr2[,model$gene]
pred.y_test<-expr_test2%*%model$coef
pred.y_test<-as_tibble(data.frame(cbind(sample=rownames(pred.y_test),pred.y_test)))
colnames(pred.y_test)<-c("cellline","EDRGS")


dir.create("output")
write_tsv(pred.y_test,"output/GDSC1_ESCC_cellline_targetgene_score.txt")
# colnames(pred.y_test)<-c("cellline","EDRGS")

####====================================================================================
####======================draw EDRGS high/low group CDK4/6 expression level=============
####====================================================================================
library(tidyverse)
ESCC_cellline<-read_tsv("cellline_NTP.txt")
GDSC1_cellline<-ESCC_cellline %>% dplyr::filter(methods=="GDSC1")


GDSC1_Palbociclib<-GDSC1_Res[GDSC1_cellline$cellline,"Palbociclib_1054"]
GDSC1_Palbociclib<-data.frame(cellline=names(GDSC1_Palbociclib),Palbociclib=GDSC1_Palbociclib)

CDK4_CDK6_GDSC1_expr<-t(ESCC_GDSC1_Expr[c("CDK4","CDK6"),])
CDK4_CDK6_GDSC1_expr2<-data.frame(cellline=rownames(CDK4_CDK6_GDSC1_expr),CDK4_CDK6_GDSC1_expr)

dt<-inner_join(GDSC1_Palbociclib,pred.y_test)
dt<-inner_join(dt,CDK4_CDK6_GDSC1_expr2)
dt<-dt[,-1]
dt$Palbociclib<-ifelse(dt$Palbociclib>median(dt$Palbociclib),0,1)

data<-dt[,c("EDRGS","CDK4","CDK6","Palbociclib")]
data<-as.data.frame(data)
data$Palbociclib<-factor(data$Palbociclib,levels = c(0,1))
data$EDRGS<-as.numeric(data$EDRGS)
data$CDK4<-as.numeric(data$CDK4)
data$CDK6<-as.numeric(data$CDK6)
library(pROC)

dir.create("output/")

pdf("output/EDRGS_CDK4_CDK6_compare_Palbociclib_ROC_in_GDSC1_add_CI.pdf",width=5,height=5)
rocobj1 <- plot.roc(data$Palbociclib, data[,1],main="", percent=TRUE, col="#FE6D73")
rocobj2 <- lines.roc(data$Palbociclib, data[,2], percent=TRUE, col="#FFCB77")
rocobj3 <- lines.roc(data$Palbociclib, data[,3], percent=TRUE, col="#7770b4")


ci1<-ci(rocobj1)
ci2<-ci(rocobj2)
ci3<-ci(rocobj3)



test2 <- roc.test(rocobj1,rocobj2,method="bootstrap")
test3 <- roc.test(rocobj1,rocobj3,method="bootstrap")


labels<-c(
  paste0(colnames(data)[1],", ","AUC=",sprintf("%0.3f",rocobj1$auc/100),"(",sprintf("%0.3f",ci1[1]/100),"-",sprintf("%0.3f",ci1[3]/100),")"),
  paste0(colnames(data)[2],", ","AUC=",sprintf("%0.3f",rocobj2$auc/100),"(",sprintf("%0.3f",ci2[1]/100),"-",sprintf("%0.3f",ci2[3]/100),")",", P-value=",sprintf("%0.3f",test2$p.value)),
  paste0(colnames(data)[3],", ","AUC=",sprintf("%0.3f",rocobj3$auc/100),"(",sprintf("%0.3f",ci3[1]/100),"-",sprintf("%0.3f",ci3[3]/100),")",", P-value=",sprintf("%0.3f",test3$p.value))
)



legend("bottomright", legend=labels, 
       col=c("#FE6D73","#FFCB77","#7770b4","#227C9D", "#17C3B2","#FFCB77","#ff9966","#fdee00","#fae7b5","#f4c2c2",
             "#57b046","#b15929","#b294c6","#98777b","#7770b4",
             "red"), 
       cex=0.5,
       lwd=2)

dev.off()


