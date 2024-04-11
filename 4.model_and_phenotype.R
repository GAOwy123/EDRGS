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
# targetgene_Score$targetgene_Score<-ifelse(targetgene_Score$targetgene_Score>median(targetgene_Score$targetgene_Score),"High","Low")
df<-inner_join(df,targetgene_Score)

# df$targetgene_Score<-factor(df$targetgene_Score,levels = c("Low","High"))
df$Age<-ifelse(df$Age < 60,"<60",">=60")
df$Age<-factor(df$Age,levels = c("<60",">=60"))

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


df$Metastatics<-ifelse(df$N_info==0,"No","Yes")
df$Metastatics<-factor(df$Metastatics,levels = c("No","Yes"))
df$Status<-gsub("1","Dead",df$Status)
df$Status<-gsub("0","Alive",df$Status)
df$Status<-factor(df$Status,levels = c("Alive","Dead"))
df<-df[,c("sample","Targetgene_Score","Status","Gender","Age","Recurrence_status","Metastatics","Grade","TNM",
          "Location","Smoking_history",
          "Drinking_history")]
#=====================================================================================
#====================   targetgene Score & phenotype     =========================
#=====================================================================================
dir.create("output/18.Targetgene_Score&phenotype")
dir.create("output/18.Targetgene_Score&phenotype/phenotype")
variants<-c("Status","Gender","Age","Recurrence_status","Metastatics","Grade","TNM",
            "Location","Smoking_history",
            "Drinking_history")

for(i in variants){
  # i<-variants[1]
  sub<-df[,c("Targetgene_Score",i)]
  sub<-na.omit(sub)
  colnames(sub)<-c("Targetgene_Score","subtype")
  
  sub$subtype<-factor(sub$subtype)
  combination<- combn(unique(unlist(sub[,"subtype"])),2)
  my_comparision<-list()
  for(j in 1:ncol(combination)){
    my_comparision<-c(my_comparision,list(as.vector(combination[,j])))
  }
  
  
  p<-ggboxplot(sub,x = "subtype", y = "Targetgene_Score",
               color = "subtype", group="subtype",
               xlab = i,ylab="Targetgene Score",title="trainset",
               palette =c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2"),
               add = "jitter")+
    stat_compare_means(comparisons = my_comparision,method="t.test",label = "p.format")+#p.signif
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf(paste0("output/18.Targetgene_Score&phenotype/phenotype/",i,"_Targetgene_Score.pdf"),width=4,height=4)
  print(p)
  dev.off()
  
  tiff(paste0("output/18.Targetgene_Score&phenotype/phenotype/",i,"_Targetgene_Score.tiff"),width=400,height=400)
  print(p)
  dev.off() 
  
}



library(ggpubr)
###==============fisher test =================
df$Targetgene_Score<-ifelse(df$Targetgene_Score>median(df$Targetgene_Score),"High","Low")
df$Targetgene_Score<-factor(df$Targetgene_Score,levels = c("Low","High"))

library(reshape2)
library(ggpubr)

dir.create("output/18.Targetgene_Score&phenotype")
for(state in variants){
  # state<-variants[1]
  sub_df<-df[,c("Targetgene_Score",state)]
  sub_df<-sub_df %>% filter(sub_df[,2] !="loss")
  sub_df<-sub_df %>% drop_na()
  
  genename<-"Targetgene_Score"
  gene_smoke<-table(filter(sub_df[,c(genename,state)],sub_df[,c(genename,state)][,2]!="loss"))
  # nrow<-length(unique(df[,genename]))
  # ncol<-length(unique(df[,state][!is.na(df[,state])]))
  nrow<-nrow(gene_smoke)
  ncol<-ncol(gene_smoke)
  df_smoke = as.data.frame(prop.table(gene_smoke,1))
  pvalue <- fisher.test(matrix(gene_smoke,ncol=ncol,nrow=nrow),simulate.p.value=TRUE) #
  
  write.table(cbind(genename,state,round(pvalue$p.value,5)),"output/18.Targetgene_Score&phenotype/pheno_fisher_test_pvalue.txt",sep="\t",append = T,quote=F,row.names = F,col.names = F)
  p_smoke =ggplot(data=df_smoke, aes(x=df_smoke[,genename], y=Freq, fill=df_smoke[,state])) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#8ECFC9","#FFBE7A","#FA7F6F","#82B0D2","#BEB8DC","#E7DAD2"))+
    # geom_text(aes(label=paste0(round(Freq*100,1),"%")))+
    geom_text(aes(label=paste0(round(Freq*100,1),"%")),
              size=5,position = position_stack(vjust = 0.5))+ #调整文本位置
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0,color="black",size=12,hjust=0.5,lineheight=1),
          axis.text.y = element_text(color="black",size=12,hjust=1,lineheight=1),
          axis.title = element_text(size=16))+
    labs(title=paste(state,paste0("P value = ",round(pvalue$p.value,3)),sep="   "),y="Frequence",x="Targetgene Score")+
    ##移除背景色和网格，加坐标轴
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),
          # legend.title=element_blank(),
          legend.position = "top",
          plot.title = element_text(hjust=0.5))+
    guides(fill = guide_legend(title = state))
  
  # p_smoke
  pdf(paste("output/18.Targetgene_Score&phenotype/fisherTest_",genename,"_",state,".pdf",sep=""), width=4, height=4.5)
  print(p_smoke)
  dev.off()
  tiff(paste("output/18.Targetgene_Score&phenotype/fisherTest_",genename,"_",state,".tiff",sep=""), width=400, height=450)
  print(p_smoke)
  dev.off()
}




