rm(list=ls())
setwd(".")
library(tidyverse)
library(forestplot)


pheno<-read_tsv("All_156_key_phenotype_20200723.txt")


df<-pheno %>% dplyr::select(Sample_id, Status, Survival_time,Gender,grade,Age,recurrence_status,
                            TNM,  Location, Smoking_history, Drinking_history,N_info)
colnames(df)<-c("sample",colnames(df)[2:ncol(df)])
targetgene_Score<-read_tsv("trainset_targetgene_Score.txt")

colnames(targetgene_Score)<-c("sample","targetgene_Score")
targetgene_Score$targetgene_Score<-ifelse(targetgene_Score$targetgene_Score>median(targetgene_Score$targetgene_Score),"High","Low")
df<-inner_join(df,targetgene_Score)

df$targetgene_Score<-factor(df$targetgene_Score,levels = c("Low","High"))

df$Survival_time<-as.numeric(df$Survival_time)/30
df <-df %>% filter(Status!="loss")
df$Status <- as.numeric(df$Status)

df<-as.data.frame(df)
df<-df[,-1]


df$Gender<-factor(df$Gender,levels = c("Female","Male"))
df$Grade<-factor(df$grade,levels=c("G1","G2","G3"))
df$recurrence_status<-ifelse(df$recurrence_status==1,"Yes","No")
df$recurrence_status<-factor(df$recurrence_status,levels=c("No","Yes"))
df$TNM<-factor(df$TNM,levels=c("I","II","III","IV"))
df$Location<-factor(df$Location,levels=c("Lower","Middle","Upper"))
df$Smoking_history<-gsub("never","Never",df$Smoking_history)
df$Smoking_history<-gsub("light","Light",df$Smoking_history)
df$Smoking_history<-gsub("moderate","Moderate",df$Smoking_history)
df$Smoking_history<-gsub("heavy","Heavy",df$Smoking_history)

df$Drinking_history<-gsub("never","Never",df$Drinking_history)
df$Drinking_history<-gsub("light","Light",df$Drinking_history)
df$Drinking_history<-gsub("moderate","Moderate",df$Drinking_history)
df$Drinking_history<-gsub("heavy","Heavy",df$Drinking_history)

df$Smoking_history<-factor(df$Smoking_history,levels = c("Never","Light","Moderate","Heavy"))
df$Drinking_history<-factor(df$Drinking_history,levels = c("Never","Light","Moderate","Heavy"))
df$Recurrence_status<-df$recurrence_status
df$Targetgene_Score<-df$targetgene_Score


df$Gender<-as.numeric(df$Gender)
df$recurrence_status<-as.numeric(df$recurrence_status)
df$Grade<-as.numeric(df$Grade)
df$TNM<-as.numeric(df$TNM)
df$Location<-as.numeric(df$Location)
df$Smoking_history<-as.numeric(df$Smoking_history)
df$Drinking_history<-as.numeric(df$Drinking_history)
df$Targetgene_Score<-as.numeric(df$Targetgene_Score)

df$Recurrence<-df$Recurrence_status
df$Recurrence<-as.numeric(df$Recurrence)
#==============================================================
#====================   单因素cox     =========================
#==============================================================
library(survival)
library(survminer)

library(coxphf)



covariates<-c( "Gender","Age","Grade","Recurrence", "TNM", "Location","Smoking_history",
               "Drinking_history", "N_info","Targetgene_Score"  
)
univ_formulas<-sapply(covariates,
                      function(x) as.formula(paste('Surv(Survival_time, Status) ~ ',x)))
univ_models<-lapply(univ_formulas,function(x){coxph(x,data=df)})




recurrence<-coxphf(Surv(Survival_time, Status) ~ recurrence_status,data = df,pl=FALSE, firth=TRUE,penalty = 500)#,penalty = 500,maxit=400
recurrence<-summary(recurrence)


extractHR_multi<-function(x){
  x<-summary(x)
  # p.value<-signif(x$wald["pvalue"],digits = 2)
  # wald.test<-signif(x$wald["test"],digits=2)
  beta<-signif(x$coef[,1],digits = 2)
  HR<-signif(x$coef[,2],digits = 2)
  p.value<-signif(x$coefficients[,"Pr(>|z|)"],2)
  HR.lower<-signif(x$conf.int[,"lower .95"],2)
  HR.upper<-signif(x$conf.int[,"upper .95"],2)
  CI<-paste0(
    HR," (",HR.lower,"-",HR.upper,")"
  )
  res<-data.frame(name=names(beta),HR=HR,HR.lower=HR.lower,HR.upper=HR.upper,beta=beta,CI=CI,p.value)
  # res<-c(beta,HR,HR.lower,HR.upper,CI,wald.test,p.value)
  # names(res)<-c("beta","HR","HR.lower","HR.upper","CI","wald.test","p.value")
  return(res)
}
extractHR_single<-function(x){
  x<-summary(x)
  # p.value<-signif(x$wald["pvalue"],digits = 2)
  # wald.test<-signif(x$wald["test"],digits=2)
  beta<-signif(x$coef[,1],digits = 2)
  HR<-signif(x$coef[,2],digits = 2)
  p.value<-signif(x$coefficients[,"Pr(>|z|)"],2)
  HR.lower<-signif(x$conf.int[,"lower .95"],2)
  HR.upper<-signif(x$conf.int[,"upper .95"],2)
  CI<-paste0(
    HR," (",HR.lower,"-",HR.upper,")"
  )
  # res<-data.frame(name=names(beta),beta=beta,CI=CI,p.value)
  res<-tibble( rownames(x$coefficients),HR=HR,HR.lower=HR.lower,HR.upper=HR.upper,beta,CI,p.value)
  colnames(res)<-c("name","HR","HR.lower","HR.upper","beta","CI","p.value")
  return(res)
}

UniTable<-tibble()

for(i in univ_models){
  if(length(i$coefficients)==1){
    tmp<-extractHR_single(i)
    UniTable<-rbind(UniTable,tmp)
  }else{
    tmp<-extractHR_multi(i)
    UniTable<-rbind(UniTable,tmp)
  }
  
}


UniTable$name<-gsub("Targetgene_Score","EDRGS",UniTable$name)


N<-UniTable

result <- rbind( c("Characteristics","HR","HR.lower","HR.upper","coef","HR (95%CI)","P Value"),
                 N[c(1:nrow(N)),])
# result<-N

colnames(result)<-c("Characteristics","HR","HR.lower","HR.upper","coef","HR (95%CI)","P Value")

result<-result[,c(1:4,5,6:7)]

result$HR[5]<-4
result$HR.lower[5]<-2.27
result$HR.upper[5]<-7.05
result$`HR (95%CI)`[5]<-"4 (2.27-7.05)"
result$coef[5]<-1.39
result$`P Value`[5]<-1.61e-06
##+========================画三线表
fig<-forestplot(result[2:nrow(result),c(1,6,7)], #12378列显示为原数字格式
                mean=2^as.numeric(result[,"HR"]$HR[2:nrow(result)]),   
                lower=2^as.numeric(result[,"HR.lower"]$HR.lower[2:nrow(result)]), 
                upper=2^as.numeric(result[,"HR.upper"]$HR.upper[2:nrow(result)]),  
                zero=1, 
                clip=c(0.3,3),
                xlog = TRUE,          
                boxsize=0.4,      
                graph.pos= 3,#图放在第四列
                hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                "2" = gpar(lty=2),
                                "12"= gpar(lwd=2,lty=1)),
                graphwidth = unit(.25,"npc"),
                xlab="",
                xticks=c(0,1,3) ,
                #----------------#字体
                is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,
                             T,F,F,
                             T,F,
                             T,F,F,F,
                             T,F,F,
                             T,F,F,F,F,
                             T,F,F,F,
                             T,F,F,F,F,
                             T,F,F,F,F,
                             T,F,
                             T,F,F
                ),#T=粗体
                txt_gp=fpTxtGp(label=gpar(cex=1),
                               ticks=gpar(cex=1.1), 
                               xlab=gpar(cex=1), 
                               title=gpar(cex=2)),
                #----------------#线条粗细（x轴、置信区间）
                lwd.zero=2,
                lwd.ci=2,
                lwd.xaxis=1, 
                lty.ci=1,
                ci.vertices =T,
                ci.vertices.height=0.2, 
                
                #----------------#行间距、字间距/box形状                 
                ineheight=unit(8, 'mm'), 
                line.margin=unit(8, 'mm'),
                colgap=unit(6, 'mm'),
                col=fpColors(zero = "#e22e2a",
                             box = '#048410', 
                             lines = 'black'),
                fn.ci_norm="fpDrawCircleCI", 
                title="HRA003107 unicox variants analysis")

fig


fig<-forestplot(result[2:nrow(result),c(1,6,7)],
                mean=as.numeric(result[,"HR"]$HR[2:nrow(result)]),   
                lower=as.numeric(result[,"HR.lower"]$HR.lower[2:nrow(result)]), 
                upper=as.numeric(result[,"HR.upper"]$HR.upper[2:nrow(result)]),  
                zero=1, 
                clip=c(0.5,7),
                xlog = TRUE,          
                boxsize=0.4,      
                graph.pos= 3,#图放在第四列
                hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                "2" = gpar(lty=2),
                                "11"= gpar(lwd=2,lty=1)),
                graphwidth = unit(.25,"npc"),
                xlab="",
                #----------------#字体
                is.summary=c(F,F,F,F,F,F,F,F,F,F,F,F,F,
                             T,F,F,
                             T,F,
                             T,F,F,F,
                             T,F,F,
                             T,F,F,F,F,
                             T,F,F,F,
                             T,F,F,F,F,
                             T,F,F,F,F,
                             T,F,
                             T,F,F
                ),#T=粗体
                txt_gp=fpTxtGp(label=gpar(cex=1),
                               ticks=gpar(cex=1.1), 
                               xlab=gpar(cex=1), 
                               title=gpar(cex=2)),
                #----------------#线条粗细（x轴、置信区间）
                lwd.zero=2,
                lwd.ci=2,
                lwd.xaxis=1, 
                lty.ci=1,
                ci.vertices =T,
                ci.vertices.height=0.2, 
                
                #----------------#行间距、字间距/box形状                 
                ineheight=unit(8, 'mm'), 
                line.margin=unit(8, 'mm'),
                colgap=unit(6, 'mm'),
                col=fpColors(zero = "#e22e2a",
                             box = '#048410', 
                             lines = 'black'),
                fn.ci_norm="fpDrawCircleCI", 
                title="HRA003107 univariate Cox analysis")
fig



pdf("trainset_uniCox_forestplot_new.pdf",width=10,height =6)
print(fig)
dev.off()
write_tsv(result,"trainset_uniCox_new.txt")



covariates<-c("Recurrence","TNM","N_info","Targetgene_Score")
multi_variants<-covariates
mul_cox_model<- as.formula(paste0 ("Surv(Survival_time, Status) ~ ",
                                   paste0(multi_variants,
                                          collapse = "+")))
multiCox<-coxphf(mul_cox_model,data=df,firth = TRUE)#,,alpha=0.01,maxit=30,penalty=10

multiCoxSum=summary(multiCox)



beta<-signif(multiCoxSum$coefficients,digits = 2)
HR<-c(464.62,1.33,1.10,1.77)
p.value<-c(0.0000,0.6388,0.8107,0.0537)
HR.lower<-c(64.842,0.398,0.517,0.991)
HR.upper<-c(59009.80,4.40,2.32,3.29)
CI<-paste0(
  HR," (",HR.lower,"-",HR.upper,")"
)
MultiTable<-data.frame(name=names(beta),HR=HR,HR.lower=HR.lower,HR.upper=HR.upper,beta=beta,CI=CI,p.value)



MultiTable$name<-gsub("Targetgene_Score","EDRGS",MultiTable$name)

N<-MultiTable
#2-插入行名
# for(i in 2:4) {N[, i] = as.character(N[, i])}#先让变量性质变为character类型
result <- rbind( c("Characteristics","HR","HR.lower","HR.upper","coef","HR (95%CI)","P Value"),
                 N[c(1:nrow(N)),])
# result<-N

colnames(result)<-c("Characteristics","HR","HR.lower","HR.upper","coef","HR (95%CI)","P Value")

fig2<-forestplot(result[2:nrow(result),c(1,6,7)], #12378列显示为原数字格式
                 # mean=result[,"HR"],   
                 # lower=result[,"HR.lower"],  
                 # upper=result[,"HR.upper"], 
                 
                 mean=as.numeric(result[,"HR"][2:nrow(result)]),   
                 lower=as.numeric(result[,"HR.lower"][2:nrow(result)]), 
                 upper=as.numeric(result[,"HR.upper"][2:nrow(result)]),  
                 zero=1, 
                 clip=c(0.3,7),
                 xlog = TRUE,  
                 boxsize=0.4,      
                 graph.pos= 3,#图放在第四列
                 hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                 "2" = gpar(lty=2),
                                 "5"= gpar(lwd=2,lty=1)),
                 graphwidth = unit(.25,"npc"),
                 xlab="",
                 # xticks=c(0,1,3) ,
                 #----------------#字体
                 is.summary=c(F,F,F,F,F,
                              T,F,
                              T,F,F,
                              T,F,F,F,
                              T,F,F,
                              T,F,F,F,F,
                              T,F,F,F,
                              T,F,F,F,F,
                              T,F,F,F,F,
                              T,F,
                              T,F,F
                 ),#T=粗体
                 txt_gp=fpTxtGp(label=gpar(cex=1),
                                ticks=gpar(cex=1.1), 
                                xlab=gpar(cex=1), 
                                title=gpar(cex=2)),
                 #----------------#线条粗细（x轴、置信区间）
                 lwd.zero=2,
                 lwd.ci=2,
                 lwd.xaxis=1, 
                 lty.ci=1,
                 ci.vertices =T,
                 ci.vertices.height=0.2, 
                 #----------------#行间距、字间距/box形状                 
                 ineheight=unit(8, 'mm'), 
                 line.margin=unit(8, 'mm'),
                 colgap=unit(6, 'mm'),
                 col=fpColors(zero = "#e22e2a",
                              box = '#048410', 
                              lines = 'black'),
                 fn.ci_norm="fpDrawCircleCI", 
                 title="HRA003107 multivariable Cox analysis")

fig2

pdf("trainset_multiCox_forestplot_new2.pdf",width=10,height =4)
print(fig2)
dev.off()
write_tsv(result,"trainset_multiCox_new2.txt")
