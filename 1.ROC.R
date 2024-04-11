rm(list=ls())
setwd(".")
library(tidyverse)
hubgene<-read_tsv("hubgene.txt") %>% unlist()

gene_expr<-read_tsv("escc.phaseII.hisat2.genesymbol.median.log2.TPM+1.txt")
mRNA<-read.table("protein_coding_gene.txt")
gene_expr<-gene_expr %>% filter(gene %in% mRNA$gene_name)
gene_expr2<-gene_expr[,-1]
gene_expr2<-as.data.frame(gene_expr2)
rownames(gene_expr2)<-gene_expr$gene
expr<-gene_expr2[hubgene,]


expr<-t(expr)
names<-rownames(expr)
expr<-as_tibble(expr)
expr$accession<-names


group<-ifelse(grepl("T",expr$accession),1,0)
pheno<-data.frame(accession=expr$accession,group=group)
data<-inner_join(pheno,expr)
samples<-data.frame(accession=data[,1])
samples$SampleType<-ifelse(grepl("T",samples$accession),"Tumor","Normal")
samples$accession<-gsub("^X","",samples$accession)
data<-data[,-1]

library(glmnet)
#data$group<-as.factor(data$group)
model1 <- glm(data$group ~ .,
              data=data, family='binomial')

pre <- predict(model1,type='response')


eRNA_score<-data.frame(samples,pre)

library(ggpubr)
eRNA_score$SampleType<-factor(eRNA_score$SampleType,levels=c("Normal","Tumor"))
p<-ggboxplot(eRNA_score,x='SampleType',y="pre",fill = "SampleType",
             palette = c("#cbdb8f",  "#eb9ea1","#d5bdd7"),main="trainset",
             ylab="eRNA Score",xlab="")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+
  stat_compare_means(method="wilcox.test",label="p.signif",label.x=1.5)





options(digits = 3)
library(pROC)
pdf("hubgene_ROC.pdf",width=6,height=6)

rocobj1 <- plot.roc(data$group, data[,2],main="", percent=TRUE, col="#1c61b6")
rocobj2 <- lines.roc(data$group, data[,3], percent=TRUE, col="#008600")
rocobj3 <- lines.roc(data$group, data[,4], percent=TRUE, col="#e9d66b")
rocobj4 <- lines.roc(data$group, data[,5], percent=TRUE, col="#ff9966")
rocobj5 <- lines.roc(data$group, data[,6], percent=TRUE, col="#fdee00")
rocobj6 <- lines.roc(data$group, data[,7], percent=TRUE, col="#fae7b5")
rocobj7 <- lines.roc(data$group, data[,8], percent=TRUE, col="#f4c2c2")
rocobj8 <- lines.roc(data$group, data[,9], percent=TRUE, col="#57b046")
rocobj9 <- lines.roc(data$group, data[,10], percent=TRUE, col="#b15929")
rocobj10 <- lines.roc(data$group, data[,11], percent=TRUE, col="#b294c6")
rocobj11 <- lines.roc(data$group, data[,12], percent=TRUE, col="#98777b")
rocobj12 <- lines.roc(data$group, data[,13], percent=TRUE, col="#7770b4")
rocobj13 <- lines.roc(data$group, pre, percent=TRUE, col="red")



ci1<-ci(rocobj1)
ci2<-ci(rocobj2)
ci3<-ci(rocobj3)
ci4<-ci(rocobj4)
ci5<-ci(rocobj5)
ci6<-ci(rocobj6)
ci7<-ci(rocobj7)
ci8<-ci(rocobj8)
ci9<-ci(rocobj9)
ci10<-ci(rocobj10)
ci11<-ci(rocobj11)
ci12<-ci(rocobj12)

ci13<-ci(rocobj13)


test1 <- roc.test(rocobj13,rocobj1,method="bootstrap")
# text(65,100,labels = paste("Combined vs ACVR1C P-value: ",round(test1$p.value,3))) # 将比较结果p值展示在图中

test2 <- roc.test(rocobj13,rocobj2,method="bootstrap")
# text(65,95,labels = paste("EDRGS vs TIDE P-value: ",round(test2$p.value,3))) # 将比较结果p值展示在图中

test3 <- roc.test(rocobj13,rocobj3,method="bootstrap")
# text(65,90,labels = paste("PDL1 vs TIDE P-value: ",round(test3$p.value,3))) # 将比较结果p值展示在图中
test4 <- roc.test(rocobj13,rocobj4,method="bootstrap")
test5 <- roc.test(rocobj13,rocobj5,method="bootstrap")
test6 <- roc.test(rocobj13,rocobj6,method="bootstrap")
test7 <- roc.test(rocobj13,rocobj7,method="bootstrap")

test8 <- roc.test(rocobj13,rocobj8,method="bootstrap")
test9 <- roc.test(rocobj13,rocobj9,method="bootstrap") 
test10 <- roc.test(rocobj13,rocobj10,method="bootstrap")
test11 <- roc.test(rocobj13,rocobj11,method="bootstrap") 
test12 <- roc.test(rocobj13,rocobj12,method="bootstrap")

labels<-c(
  paste0(colnames(data)[2],", ","AUC=",sprintf("%0.3f",rocobj1$auc/100),"(",sprintf("%0.3f",ci1[1]/100),"-",sprintf("%0.3f",ci1[3]/100),")",", P-value=",signif(test1$p.value,3)/100),
  paste0(colnames(data)[3],", ","AUC=",sprintf("%0.3f",rocobj2$auc/100),"(",sprintf("%0.3f",ci2[1]/100),"-",sprintf("%0.3f",ci2[3]/100),")",", P-value=",signif(test2$p.value,3)/100),
  paste0(colnames(data)[4],", ","AUC=",sprintf("%0.3f",rocobj3$auc/100),"(",sprintf("%0.3f",ci3[1]/100),"-",sprintf("%0.3f",ci3[3]/100),")",", P-value=",signif(test3$p.value,3)/100),
  paste0(colnames(data)[5],", ","AUC=",sprintf("%0.3f",rocobj4$auc/100),"(",sprintf("%0.3f",ci4[1]/100),"-",sprintf("%0.3f",ci4[3]/100),")",", P-value=",signif(test4$p.value,3)/100),
  paste0(colnames(data)[6],", ","AUC=",sprintf("%0.3f",rocobj5$auc/100),"(",sprintf("%0.3f",ci5[1]/100),"-",sprintf("%0.3f",ci5[3]/100),")",", P-value=",signif(test5$p.value,3)/100),
  paste0(colnames(data)[7],", ","AUC=",sprintf("%0.3f",rocobj6$auc/100),"(",sprintf("%0.3f",ci6[1]/100),"-",sprintf("%0.3f",ci6[3]/100),")",", P-value=",signif(test6$p.value,3)/100),
  paste0(colnames(data)[8],", ","AUC=",sprintf("%0.3f",rocobj7$auc/100),"(",sprintf("%0.3f",ci7[1]/100),"-",sprintf("%0.3f",ci7[3]/100),")",", P-value=",signif(test7$p.value,3)/100),
  paste0(colnames(data)[9],", ","AUC=",sprintf("%0.3f",rocobj8$auc/100),"(",sprintf("%0.3f",ci8[1]/100),"-",sprintf("%0.3f",ci8[3]/100),")",", P-value=",signif(test8$p.value,3)/100),
  paste0(colnames(data)[10],", ","AUC=",sprintf("%0.3f",rocobj9$auc/100),"(",sprintf("%0.3f",ci9[1]/100),"-",sprintf("%0.3f",ci9[3]/100),")",", P-value=",signif(test9$p.value,3)/100),
  paste0(colnames(data)[11],", ","AUC=",sprintf("%0.3f",rocobj10$auc/100),"(",sprintf("%0.3f",ci10[1]/100),"-",sprintf("%0.3f",ci10[3]/100),")",", P-value=",signif(test10$p.value,3)/100),
  paste0(colnames(data)[12],", ","AUC=",sprintf("%0.3f",rocobj11$auc/100),"(",sprintf("%0.3f",ci11[1]/100),"-",sprintf("%0.3f",ci11[3]/100),")",", P-value=",signif(test11$p.value,3)/100),
  paste0(colnames(data)[13],", ","AUC=",sprintf("%0.3f",rocobj12$auc/100),"(",sprintf("%0.3f",ci12[1]/100),"-",sprintf("%0.3f",ci12[3]/100),")",", P-value=",signif(test12$p.value,3)/100),
  # paste0("TMEM158, ","AUC=",signif(rocobj7$auc,3)/100),
  # paste0("GSN, ","AUC=",signif(rocobj8$auc,3)/100),
  paste0("Combined, ","AUC=",sprintf("%0.3f",rocobj13$auc/100),"(",sprintf("%0.3f",ci13[1]/100),"-",sprintf("%0.3f",ci13[3]/100),")")
)



legend("bottomright", legend=labels, 
       col=c("#1c61b6", "#008600","#e9d66b","#ff9966","#fdee00","#fae7b5","#f4c2c2",
             "#57b046","#b15929","#b294c6","#98777b","#7770b4",
             "red"), 
       lwd=2,cex=0.75)

dev.off()





##计算ROC曲线
roc_obj <- roc(data$group, pre) 

##计算AUC的置信区间
ci<-ci(roc_obj)
cat("AUC:", auc(roc_obj),"\n")
cat("95% Confidence Interval: ", ci[1],"-",ci[3],"\n")

# 95% Confidence Interval:  0.9976318 - 1 

library(caret)
prediction_res<-ifelse(pre>0.5,1,0)
data$group<-factor(data$group,levels = c(0,1))
prediction_res<-factor(prediction_res,levels = c(0,1))
conf_mat<-confusionMatrix(prediction_res,data$group,positive = "1")

sensitivity<-conf_mat$byClass["Sensitivity"]
specificity<-conf_mat$byClass["Specificity"]






