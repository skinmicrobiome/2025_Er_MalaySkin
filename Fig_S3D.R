
library(readxl)
library(knitr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(microbiome)
library(stringr)
library(phyloseq)
library(vegan)
library(pacman)

pacman::p_load( clusterProfiler,MicrobiomeStat,tidyverse, dplyr, ggplot2,ggsci,vegan,reshape2,phyloseq,ggpubr,readr,tidyverse)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


data5=read.csv("kegg_cory_2805.csv")


levels(data5$oauc) = c("Orang Asli", "Urban")
t.test2 <- lapply(data5[,c(1:755)], function(x) {
  summary(lm(x ~ data5$oauc, var.equal=FALSE))$coefficients[2,4]
})
p_values <- unlist(t.test2)
w=data.frame(p_values)


t.test3 <- lapply(data5[,c(1:755)], function(x) {
  summary(lm(x ~ data5$oauc, var.equal=FALSE))$coefficients[2,1]
})
estimates <- unlist(t.test3)

est <- data.frame(estimates)


w=cbind(w,est)

w2=filter(w,p_values<0.05)
w2$oauc <- ifelse(w2$estimates > 0, "OA", "UC")

w2$pathway=rownames(w2)

c3=filter(w2, oauc =="OA")
c4=filter(w2, oauc =="UC")


n <- data.frame(cohort = "OA", pathway = paste(c3$pathway, collapse = ","))
n2 <- data.frame(cohort = "UC", pathway = paste(c4$pathway, collapse = ","))


eggnog.df=rbind(n,n2)

uniq_kegg <- unique(unlist(lapply(eggnog.df$pathway, function(x){unlist(str_split(x, ","))})))

eggnog.df=n

kk_all <- data.frame()
for (i in 1:dim(eggnog.df)[1]){
    l <- unique(unlist(lapply(eggnog.df$pathway, function(x){unlist(str_split(x, ","))})))
    kk <- enrichKEGG(unique(l), organism = 'ko', keyType = "kegg", pvalueCutoff = 0.01, pAdjustMethod="BH")
    kk <- as.data.frame(kk)
    kk$genome <- eggnog.df$cohort
    kk_all <- rbind(kk_all, kk)
}

eggnog.df2=n2

kk_all2 <- data.frame()
for (i in 1:dim(eggnog.df2)[1]){
    l <- unique(unlist(lapply(eggnog.df2$pathway, function(x){unlist(str_split(x, ","))})))
    kk <- enrichKEGG(unique(l), organism = 'ko', keyType = "kegg", pvalueCutoff = 0.01, pAdjustMethod="BH")
    kk <- as.data.frame(kk)
    kk$genome <- eggnog.df2$cohort
    kk_all2 <- rbind(kk_all2, kk)
}


kk_all=rbind(kk_all,kk_all2)


kk_all=filter(kk_all, !Description =="Vancomycin resistance")

p2=ggplot(kk_all,
       aes(x=-log10(p.adjust), y=Description, color=genome, size=Count)) +
geom_point(alpha=0.7) +scale_color_manual(values=c("seagreen","sienna3"))+ theme_classic()+labs(x="-log10(q-value)",y="KEGG Pathway",size="Count")+
  theme( axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x =  element_text(size=8),
         plot.title=element_text(size=8),legend.position="bottom")+
  guides(color = guide_legend(nrow=2,override.aes = list(size=5)),size = guide_legend(nrow=2))


ggsave("Fig_S3D.svg",p2+theme(legend.position="none"),dpi=320,height=3,width=4)


