
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

pacman::p_load(tidyverse, dplyr, ggplot2,ggsci,vegan,reshape2,phyloseq,ggpubr,readr,tidyverse)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")





ps=readRDS("phyloseq_0505.rda")
tax=data.frame(tax_table(ps))
tax$sp2=paste(tax$Genus,tax$Species)
tax2=tax_table(tax)
taxa_names(tax2)=rownames(tax)
ps=phyloseq(tax2, otu_table(ps),sample_data(ps))

ps=microbiome::transform(ps, "Z")
colnames(tax_table(ps))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species","sp2")


ps2=subset_taxa(ps, sp2 %in% c("Moraxella sp.","Micrococcus luteus","Corynebacterium tuberculostearicum",
"Cutibacterium acnes","Staphylococcus epidermidis","Dietzia sp.","Kocuria palustris"
))




df3=psmelt(ps2)


colors = c("#FFBC00","#FF8000","#b1ffb1",
                          "#27ff27","#00b100","#25A18E","#004e00")



df4=filter(df3, site_specific %in% c("Fh"))
df4=filter(df4, sp2 %in% c("Cutibacterium acnes"))
df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))


a <- ggboxplot(df4, x ="village" , y = "Abundance",
               color = "village",               )+theme_classic()+
  labs(x="Forehead",y="log(1+Abundance)",
       color="Status")+  
scale_color_manual(values=colors)+
facet_grid(~ sp2)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size=8,face="italic"),
         legend.position="none")+ylim(-1.6,1.6)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)



df4=filter(df3, site_specific %in% c("Vf"))
df4=filter(df4, sp2 %in% c("Cutibacterium acnes"))
df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))


b <- ggboxplot(df4, x ="village" , y = "Abundance",
               color = "village",               )+theme_classic()+
  labs(x="Volar forearm",y="log(1+Abundance)",
       color="Status")+  
scale_color_manual(values=colors)+
facet_grid(~ sp2)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
         axis.title.x =element_blank(),
         axis.title.y = element_blank(),
         strip.text.x = element_text(size=8,face="italic"),
legend.position="none")+ylim(-1.6,1.6)



df4=filter(df3, site_specific %in% c("Tw"))
df4=filter(df4, sp2 %in% c("Corynebacterium tuberculostearicum","Staphylococcus epidermidis"))
df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

c <- ggboxplot(df4, x ="village" , y = "Abundance",
               color = "village",              )+theme_classic()+
  labs(x="Toe web",y="log(1+Abundance)",
       color="Status")+  
scale_color_manual(values=colors)+
facet_grid(~ sp2)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         strip.text.x = element_text(size=8,face="italic"),
         )+ylim(-3.1,3.1)




q=ggarrange(a,b,c+theme(legend.position="none"),nrow=1,ncol=3,widths=c(0.5,0.4,1))



df4=filter(df3, site_specific %in% c("Fh"))
df4=filter(df4, sp2 %in% c("Micrococcus luteus"))
df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

a <- ggboxplot(df4, x ="village" , y = "Abundance",
               color = "village",               )+theme_classic()+
  labs(x="Forehead",y="log(1+Abundance)",
       color="Status")+  
scale_color_manual(values=colors)+
facet_grid(~ sp2)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size=8,face="italic"),
         legend.position="none")+ylim(-3.5,2)



df4=filter(df3, site_specific %in% c("Vf"))
df4=filter(df4, sp2 %in% c("Dietzia sp."))
df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

b <- ggboxplot(df4, x ="village" , y = "Abundance",
               color = "village",               )+theme_classic()+
  labs(x="Volar forearm",y="log(1+Abundance)",
       color="Status")+  
scale_color_manual(values=colors)+
facet_grid(~ sp2)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
         axis.title.x =element_text(size=8),
         axis.title.y = element_blank(),
         strip.text.x = element_text(size=8,face="italic"),
legend.position="none")+ylim(-2,2)



df4=filter(df3, site_specific %in% c("Tw"))
df4=filter(df4, sp2 %in% c("Kocuria palustris","Moraxella sp."))
df4$village=factor(df4$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

c <- ggboxplot(df4, x ="village" , y = "Abundance",
               color = "village",              )+theme_classic()+
  labs(x="Toe web",y="log(1+Abundance)",
       color="Status")+  
scale_color_manual(values=colors)+
facet_grid(~ sp2)+
 stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Kuala Lumpur",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
axis.ticks.y = element_blank(),
         axis.title.x = element_text(size=8),
         axis.title.y = element_blank(),
         strip.text.x = element_text(size=8,face="italic"),
         )+ylim(-3.5,3)




q2=ggarrange(a,b,c+theme(legend.position="none"),nrow=1,ncol=3,widths=c(0.5,0.4,1))


q3=ggarrange(q,q2,nrow=2,ncol=1,heights=c(0.9,1))

c=c+theme(legend.position="bottom")+
  guides(color = guide_legend(override.aes = list(size=5),nrow=1))


q4=ggarrange(q3,get_legend(c),nrow=2,ncol=1,heights=c(1,0.1))

ggsave("Fig_S4C-D.svg",q4,height=4.1,width=8)





