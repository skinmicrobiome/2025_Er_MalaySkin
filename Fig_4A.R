
library(readxl)
library(knitr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(XML)
library(RcppParallel)
library(stringr)
library(phyloseq)
library(vegan)
library(rstatix)
library(ggpubr)



setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")




ps=readRDS("phyloseq_0505.rda")

meta=data.frame(sample_data(ps))
df=data.frame(library_id=meta$library_id, subject_id=meta$subject_id,site=meta$site_specific,
              Disease=meta$Disease, AU=meta$a_u, village=meta$village,gender=meta$Gender,age=meta$Group,BMI=meta$Category,pets=meta$Pets_1)


otu=data.frame((otu_table(ps)))
shannon=diversity(t(otu), index="shannon")
shannon=data.frame(shannon)
shannon$library_id=rownames(shannon)
final=dplyr::left_join(shannon,df)

write.csv(final,"alpha_180624.csv")

df4=read.csv("alpha_180624.csv")

df4=filter(df4, site !="Ctrl")

df4$site=factor(df4$site,levels=c("Fh","Vf","Tw"))
levels(df4$site)=c("Forehead","Volar forearm","Toe web")


df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

comparison=list(c("Kuala Lumpur","Village C"),c("Kuala Lumpur","Village AB"),
                c("Kuala Lumpur","Village G"),c("Kuala Lumpur","Village D"),
                c("Kuala Lumpur","Village F"),c("Washington DC","Kuala Lumpur"))

colors = c("#FFBC00","#FF8000","#b1ffb1",
                          "#27ff27","#00b100","#25A18E","#004e00")

c <- ggboxplot(df4, x ="village" , y = "shannon",
              color = "village"
               ,add="jitter")+theme_classic()+
  labs(x=" ",y="Shannon index",
       color="Location")+ scale_color_manual(values=colors)+ #stat_compare_means(aes(label = ..p.signif..),label.x=1.5,label.y=13,size=3)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)+facet_wrap(~site)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_blank(),
         plot.title=element_text(size=8,face="bold"))

ggsave("Fig_4A.svg",c+theme(legend.position="bottom"),dpi=320,height=3,width=7.7)
