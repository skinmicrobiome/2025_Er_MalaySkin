
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




tax=data.frame(tax_table(ps))
tax2=tax%>% select(1:7)
tax3=tax_table(tax2)
taxa_names(tax3)=rownames(tax2)
tax_table(ps)=tax3



colnames(tax_table(ps))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

bracken.phy=subset_samples(ps,!village %in% c("Washington DC","Kuala Lumpur"))

bracken.phy=subset_taxa(bracken.phy, Species == "concentricum")


bracken.phy=subset_samples(bracken.phy, site_specific != "Ctrl")
bracken.phy <- subset_taxa(bracken.phy, !is.na(Phylum))


df2=psmelt(bracken.phy)


df2$site_specific=factor(df2$site_specific,levels=c("Fh","Vf","Tw"))
levels(df2$site_specific)=c("Forehead","Volar forearm","Toeweb")



df2$Disease=recode(df2$Disease, "Negative"="Healthy control")
df2$status=factor(df2$a_u,levels=c("U","A"))
levels(df2$status)=c("Nonlesional skin","Lesional skin")

df2=filter(df2, Disease %in% c("Healthy control","Tinea imbricata"))

df2$group_status <- ifelse(df2$Disease == "Healthy control", "Healthy control", ifelse(df2$Disease == "Tinea imbricata", df2$status, NA))


df2$group_status=recode(df2$group_status,"2"="Lesional skin","1"="Nonlesional skin")

df2$group_status=factor(df2$group_status,levels=c("Healthy control","Nonlesional skin","Lesional skin"))




df4=filter(df2,Species =="concentricum")


comparison=list(c("Healthy control","Nonlesional skin"),c("Nonlesional skin","Lesional skin"),
                c("Healthy control","Lesional skin"))


b <- ggboxplot(df4, x ="group_status" , y = "Abundance",
               color = "group_status"
)+theme_classic()+ geom_jitter(aes(color=group_status))+
  labs(x="Status",y="Relative abundance of Trichophyton concentricum",color="Status")+scale_color_manual(values=c("royalblue","salmon","violetred3"))+
  stat_compare_means(aes(label = ..p.signif..),comparisons = comparison)+
scale_y_continuous(labels = scales::percent)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x = element_blank(),
         plot.title=element_text(size=8),legend.position="bottom")+facet_wrap(~site_specific)


ggsave("Fig_6C.svg",b,height=3,width=7)



