rm(list=ls())

setwd(".")

library(tidyverse)
library(glmnet) 
library(survival)

hubgene<-read_tsv("hubgene.txt") %>% unlist()

DEG<-hubgene

gene_expr<-read_tsv("escc.phaseII.hisat2.genesymbol.median.log2.TPM+1.txt")
mRNA<-read.table("protein_coding_gene.txt")
gene_expr<-gene_expr %>% filter(gene %in% mRNA$gene_name)
gene_expr2<-gene_expr[,-1]
gene_expr2<-as.data.frame(gene_expr2)
rownames(gene_expr2)<-gene_expr$gene
expr<-gene_expr2[DEG,]

expr<-expr[DEG,grepl("T",colnames(expr))]


library(tidyverse)
library(glmnet) ##Lasso回归
library(rms)  ## 画列线图；
library(VIM) ## 包中aggr()函数，判断数据缺失情况
library(survival) ##  生存分析包


data<-t(expr)
x<-as.matrix(data) ## row is sample, col is gene
x<-apply(x,2,as.numeric)
rownames(x)<-colnames(expr)


pheno<-read_tsv("All_156_key_phenotype_20200723.txt")
df<-pheno[,c("Status","Survival_time")]
df<-as.data.frame(df)
rownames(df)<-pheno$Sample_id
df <-df %>% filter(Status !="loss")

x<-x[rownames(df),]
# if(!all(rownames(df)==rownames(x))){
#   df<-df[rownames(x),]
# }
all(rownames(df)==rownames(x))

y<-df
y$Status<-as.double(y$Status)
y$Survival_time<-as.double(y$Survival_time)
y<-as.matrix(survival::Surv(y$Survival_time,y$Status))

##################################################################
expr<-as_tibble(cbind(sample=rownames(x),x))
df<-pheno %>% dplyr::select(Sample_id,Status,Survival_time)
colnames(df)<-c("sample","Status","Survival_time")
clinical<-inner_join(df,expr)
clinical<-clinical[,-1]
clinical[,3:ncol(clinical)]<-apply(clinical[,3:ncol(clinical)],2,as.numeric)
clinical$Status<-as.numeric(clinical$Status)
#####利用筛选出的candidates 建立多因素cox回归模型
sig_gene_multi_cox<-DEG
sig_gene_multi_cox<-ifelse(grepl("^\\d",sig_gene_multi_cox),paste0("ENSR",sig_gene_multi_cox),sig_gene_multi_cox)
colnames(clinical)<-ifelse(grepl("^\\d",colnames(clinical)),paste0("ENSR",colnames(clinical)),colnames(clinical))
sig_gene_multi_cox<-gsub(":","",sig_gene_multi_cox)
sig_gene_multi_cox<-gsub("-","",sig_gene_multi_cox)
colnames(clinical)<-gsub(":","",colnames(clinical))
colnames(clinical)<-gsub("-","",colnames(clinical))


formula_for_multivariate<-as.formula(paste0('Surv(Survival_time,Status)~',paste(sig_gene_multi_cox,sep = "",collapse = "+")))
multi_variate_cox<-coxph(formula_for_multivariate,data=clinical)
##check if variances are supported by PH hypothesis
ph_hypo_multi<-cox.zph(multi_variate_cox)
ph_hypo_table<-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]

# plot(ph_hypo_multi)
dir.create("output/10.PH_test")
setwd("output/10.PH_test")
pdf("cox_time_test.pdf",width=13,height = 10)
ggcoxzph(ph_hypo_multi)
dev.off()

dir.create("PH_test")
for(g in sig_gene_multi_cox){
  # g<-sig_gene_multi_cox[1]
  formula_for_multivariate<-as.formula(paste0('Surv(Survival_time,Status)~',paste(g,sep = "",collapse = "+")))
  multi_variate_cox<-coxph(formula_for_multivariate,data=clinical)
  ##check if variances are supported by PH hypothesis
  ph_hypo_multi<-cox.zph(multi_variate_cox)

}


## remove variances not supported by ph hypothesis and perform the 2nd regression
formula_for_multivariate<-as.formula(paste0('Surv(Survival_time,Status)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = "",collapse = "+")))
multi_variate_cox2<-coxph(formula_for_multivariate,data=clinical)


coef<-multi_variate_cox2$coefficients
final_gene<-names(coef)



model<-data.frame(final_gene,coef)
colnames(model)<-c("eRNA","coef")
write_tsv(model,"Model.txt")



data<-model %>% dplyr::arrange(desc(coef))
data$eRNA<-factor(data$eRNA,levels = data$eRNA)
pdf("coefficient.pdf",width=5,height = 5)
ggplot(data=data, mapping=aes(x=eRNA,y=coef))+
  geom_bar(stat="identity",fill=ifelse(data$coef>0,'pink','lightblue'))+
  coord_flip()+xlab("")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

dev.off()

tiff("coefficient.tiff",width = 500,height = 500)
ggplot(data=data, mapping=aes(x=eRNA,y=coef))+
  geom_bar(stat="identity",fill=ifelse(data$coef>0,'pink','lightblue'))+
  coord_flip()+xlab("")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

dev.off()


