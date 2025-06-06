
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




df5=filter(df4, !village %in% c("Washington DC","Kuala Lumpur"))

df5$Disease=recode(df5$Disease, "Negative"="Healthy control")


df6=filter(df5, Disease %in% c("Healthy control", "Tinea imbricata"))
df6$Disease=factor(df6$Disease, levels=c("Healthy control","Tinea imbricata"))

df6$group_status <- ifelse(df6$Disease == "Healthy control", "Healthy control", ifelse(df6$Disease == "Tinea imbricata", df6$AU,NA))


df6$group_status=recode(df6$group_status,"A"="Lesional skin","U"="Nonlesional skin")


df6$group_status=factor(df6$group_status,levels=c("Healthy control","Nonlesional skin","Lesional skin"))


comparison=list(c("Healthy control","Nonlesional skin"),c("Nonlesional skin","Lesional skin"),
                c("Healthy control","Lesional skin"))

df6$blank=""
q2 <- ggboxplot(df6, x ="group_status" , y = "shannon",
               color = "group_status"            )+theme_classic()+
  labs(x="Status",y="Shannon index",
       color="Status")+  
scale_color_manual(values=c("royalblue","salmon","violetred3"))+
facet_wrap(~ site,nrow=3)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Healthy control",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x =element_blank(),
         axis.title.y = element_text(size=8),
         #strip.text.x = element_blank(),
         plot.title=element_text(size=8))+coord_flip()


ggsave("Fig_S6C.svg",q2+theme(legend.position="none",axis.title.y=element_blank(),axis.text.y=element_blank()),dpi=320,height=4,width=3)

