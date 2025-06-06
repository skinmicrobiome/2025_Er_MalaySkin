
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


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


map_data=read.csv("map_data_coded.csv")

map_data=filter(map_data, site !="Ctrl")




map_data$site=factor(map_data$site,levels=c("Fh","Vf","Tw"))
levels(map_data$site)=c("Forehead","Volar forearm","Toeweb")


map_data$village=factor(map_data$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))


map_data$group=factor(map_data$group,levels=c("SMGC","SMGC+Malaysia MAGs"))

p <- ggboxplot(map_data, x = "village", y = "mapped",
               color = "group", palette =c("sienna3", "seagreen"))+theme_classic()+ 
  labs(x="Location",y="Mapping rates/ %", color="Database")+
  stat_compare_means(aes(group = group),angle=90,label = "p.signif",size=5)+
  facet_wrap(~site)+ #geom_hline(yintercept=50,color="darkblue",linewidth=1,linetype="dashed")+
  theme(axis.text.x = element_text(angle=90,hjust=0.5,vjust=0.5,size=8,face="bold"),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        axis.text.y = element_text(size=8,face="bold"),
        strip.text.x =element_text(size=8,face="bold"), 
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom",strip.background = element_rect(fill="white"))+ guides(colour = guide_legend(override.aes = list(size=5)))

ggsave("Fig_S2B.svg",p,height = 4.3,width=7.5)

