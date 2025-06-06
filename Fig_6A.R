
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


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


ps=readRDS("phyloseq_0505.rda")

ps=subset_samples(ps,!village %in% c("Kuala Lumpur","Washington DC"))

ps=subset_samples(ps, site_specific=="Fh")
w=data.frame(sample_data(ps))

n2=data.frame(table(w$village,w$Disease))
n2=filter(n2, !Var2 %in% c("Tinea","Others"))

colnames(n2)=c("village","Disease","Freq")


n2$village=factor(n2$village,levels=c("Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))


n2$Disease=factor(n2$Disease,
                levels=c("Negative","Tinea imbricata","Tinea versicolor","Scabies"
                         ))


p=ggplot(data=n2, aes(x=Freq, y=village, fill=Disease)) +
  geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c("royalblue","violetred3","thistle3","ivory3"))+
  theme_classic()+
 scale_x_continuous(labels = scales::percent) +
  theme( axis.text.x = element_text(size=8,face="bold",angle=90,vjust=0.5,hjust=0.5),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_blank(),
         plot.title=element_text(size=8,face="bold"))



ggsave("Fig_6A.svg",d+theme(legend.position="right"),dpi=320,height=2,width=7)




