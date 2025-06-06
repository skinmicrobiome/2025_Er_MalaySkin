library(pacman)

pacman::p_load(microViz,tidyverse, dplyr, microbiome,ggplot2,ggsci,vegan,phyloseq,ggpubr,readr,tidyverse)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db/input_file_r_code_coded")



ps=readRDS("phyloseq_0505.rda")

ps=subset_taxa(ps, Kingdom =="Bacteria")


tax=data.frame(tax_table(ps))
tax$sp2=paste(tax$Genus, tax$Species)
tax2=tax_table(tax)
taxa_names(tax2)=rownames(tax)
tax_table(ps)=tax2
colnames(tax_table(ps))=c("Kingdom","Phylum","Class","Order","Family","Genus",
"Species","sp2")


ps1=subset_samples(ps,site_specific != "Ctrl")
ps1=subset_samples(ps1, !village %in% c("Washington DC","Kuala Lumpur"))


ps1=subset_samples(ps1, Disease %in% c("Negative","Tinea imbricata"))




ibd <- tax_fix(ps1)

sample_data(ibd)$Disease=recode(sample_data(ibd)$Disease,"Negative"="Healthy control")



sample_data(ibd)$group_status_ti <- ifelse(sample_data(ibd)$Disease == "Healthy control", "Healthy control", ifelse(sample_data(ibd)$Disease == "Tinea imbricata", sample_data(ibd)$a_u, NA))
sample_data(ibd)$group_status_ti=recode(sample_data(ibd)$group_status_ti,"A"="Lesional skin","U"="Nonlesional skin")
sample_data(ibd)$group_status_ti=factor(sample_data(ibd)$group_status_ti,levels=c("Healthy control","Nonlesional skin","Lesional skin"))


ibd2=subset_samples(ibd, site_specific =="Fh")

p1=ibd2 %>%
  tax_transform("normalize", rank = "sp2") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "group_status_ti",plot_taxa = 1:8, size = 1, ) +
  scale_colour_manual(values=c("royalblue","salmon2","violetred3"))+theme_classic()+labs(title="Toeweb",color="Status")+
  theme(axis.text.x = element_text(size=8,face="bold"),legend.position="bottom",
        axis.text.y = element_text(size=8,face="bold"),plot.title=element_text(size=8,face="bold"),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        strip.text.x = element_text(size=8, face="bold"),
        strip.text.y = element_text(size=8, face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

#vf

ibd2=subset_samples(ibd, site_specific =="Vf")

p2=ibd2 %>%
  tax_transform("normalize", rank = "sp2") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "group_status_ti",plot_taxa = 1:3, size =1, ) +
  scale_colour_manual(values=c("royalblue","salmon2","violetred3"))+theme_classic()+labs(title="Toeweb",color="Status")+
  theme(axis.text.x = element_text(size=8,face="bold"),legend.position="bottom",
        axis.text.y = element_text(size=8,face="bold"),plot.title=element_text(size=8,face="bold"),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        strip.text.x = element_text(size=8, face="bold"),
        strip.text.y = element_text(size=8, face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))






#tw


ibd2=subset_samples(ibd, site_specific =="Tw")
#ibd3 = aggregate_taxa(ibd2, "Species")


p3=ibd2 %>%
  tax_transform("normalize", rank = "sp2") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "group_status_ti",plot_taxa = 1:4, size =1, ) +
  scale_colour_manual(values=c("royalblue","salmon2","violetred3"))+theme_classic()+labs(title="Toeweb",color="Status")+
  theme(axis.text.x = element_text(size=8,face="bold"),legend.position="bottom",
        axis.text.y = element_text(size=8,face="bold"),plot.title=element_text(size=8,face="bold"),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        strip.text.x = element_text(size=8, face="bold"),
        strip.text.y = element_text(size=8, face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))


ggsave("Fig_S6D_1.svg",p1+theme(legend.position="none"),height=3,width=2.3)

ggsave("Fig_S6D_2.svg",p2+theme(legend.position="none"),height=3,width=2.3)
ggsave("Fig_S6D_3.svg",p3+theme(legend.position="none"),height=3,width=2.3)






