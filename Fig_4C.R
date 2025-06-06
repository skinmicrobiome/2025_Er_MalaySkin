
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

library(microViz)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


ps=readRDS("phyloseq_0505.rda")

ps=subset_taxa(ps, Kingdom =="Bacteria")
#ps=transform(ps, "compositional")

sample_data(ps)$village=factor(sample_data(ps)$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

tax=data.frame(tax_table(ps))
tax$sp2=paste(tax$Genus,tax$Species)
tax2=tax_table(tax)
taxa_names(tax2)=rownames(tax)
tax_table(ps)=tax2
colnames(tax_table(ps))=c("Kingdon","Phylum","Class","Order","Family",
                          "Genus","Species","sp2")
ps=aggregate_taxa(ps,"sp2")


ibd <- tax_fix(ps)



#play around the transform with either “total”, “max”, “frequency”, “normalize”, “range”, “rank”, “rrank”, “standardize”, “pa”, “chi.square”, “hellinger”, “log”, “clr”, “rclr”, “alr”
#adjust color and shape variable accordingly

ibd1=subset_samples(ibd, site_specific =="Fh")

colors = c("#FFBC00","#FF8000","#b1ffb1",
                          "#27ff27","#00b100","#25A18E","#004e00")


p1=ibd1 %>%
  tax_transform("normalize", rank = "sp2") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "village", plot_taxa = 1:3, size = 2, tax_lab_style = tax_lab_style(
    max_angle = 90, size = 2, fontface = "bold.italic"
  )) +  scale_colour_manual(values=colors) +theme_classic()+
  theme(axis.text.y =  element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.title.y =  element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.position = "bottom",plot.title = element_text(size=8))+ guides(colour = guide_legend(override.aes = list(size=5)))

ggsave("Fig_4C_1.svg",p1+theme(legend.position="none"),height=2.3,width=2.3)

ibd2=subset_samples(ibd, site_specific =="Vf")


p2=ibd2 %>%
  tax_transform("normalize", rank = "sp2") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "village", plot_taxa = 1:3, size = 2, tax_lab_style = tax_lab_style(
    max_angle = 90, size = 3, fontface = "bold.italic"
  )) +ylim(-1,1.5)+xlim(-1,1.25)+
  scale_colour_manual(values=colors)+theme_classic()+
  theme(axis.text.y =  element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.title.y =  element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.position = "none",plot.title = element_text(size=8))+ guides(colour = guide_legend(override.aes = list(size=5)))


ggsave("Fig_4C_2.svg",p2+theme(legend.position="none"),height=2.3,width=2.3)




ibd3=subset_samples(ibd, site_specific =="Tw")


p3=ibd3 %>%
  tax_transform("normalize", rank = "sp2") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "village", plot_taxa = 1:5, size = 2, tax_lab_style = tax_lab_style(
    max_angle = 90, size = 3, fontface = "bold.italic"
  )) +labs(color="Village")+ylim(-1,1)+xlim(-1,1)+
  scale_colour_manual(values=colors) +
  coord_fixed(ratio = 1, clip = "off")+ theme_classic()+
  theme(axis.text.x = element_text(size=8),legend.position="none",
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.text.y = element_text(size = 8, face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

ggsave("Fig_4C_2.svg",p3+theme(legend.position="none"),height=2.3,width=2.3)
