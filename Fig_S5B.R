
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


df4=read.csv("alpha_180624.csv")

df4=filter(df4, site !="Ctrl")

df4$site=factor(df4$site,levels=c("Fh","Vf","Tw"))
levels(df4$site)=c("Forehead","Volar forearm","Toe web")


df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))



comparison=list(c("Forehead","Volar forearm"),c("Volar forearm","Toe web"),
                c("Forehead","Toe web"))


c <- ggboxplot(df4, x ="site" , y = "shannon",
               color = "site",
               add = "jitter")+theme_classic()+
  labs(x="Body site",y="Shannon index",
       color="Body site")+  #stat_compare_means(aes(label = ..p.signif..),label.x=1.5,label.y=10,size=3)
scale_color_manual(values=c("seagreen",
"orange3","purple3"))+stat_compare_means(label.x=1.4,label.y=10,label="p.signif")+
  stat_compare_means(aes(label = ..p.signif..),comparisons = comparison,size=3)+facet_wrap(~village,nrow=1)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x = element_blank(),legend.position="none",
         plot.title=element_text(size=8))


ggsave("Fig_S5B.svg",c+theme(axis.title.x=element_blank()),dpi=320,height=2.5,width=7.7)

