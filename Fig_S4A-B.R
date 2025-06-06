
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
library(DataExplorer)

pacman::p_load(DataExplorer ,tidyverse, dplyr, ggplot2,ggsci,vegan,reshape2,phyloseq,ggpubr,readr,tidyverse)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")



ps=readRDS("phyloseq_0505.rda")
tax=data.frame(tax_table(ps))
tax$sp2=paste(tax$Genus,tax$Species)
tax2=tax_table(tax)
taxa_names(tax2)=rownames(tax)
ps=phyloseq(tax2, otu_table(ps),sample_data(ps))
colnames(tax_table(ps))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species","sp2")


ps2=subset_taxa(ps, sp2 %in% c("Moraxella sp.","Micrococcus luteus","Corynebacterium tuberculostearicum",
"Cutibacterium acnes","Staphylococcus epidermidis","Dietzia sp.","Kocuria palustris"
))

ps2=microbiome::transform(ps2,"log")

tax=data.frame(tax_table(ps2))
otu=data.frame(t(otu_table(ps2)))
colnames(otu)=tax$sp2

sample_data(ps)$c_acnes=otu$`Cutibacterium acnes`
sample_data(ps)$s_epi=otu$`Staphylococcus epidermidis`
sample_data(ps)$c_tuber=otu$`Corynebacterium tuberculostearicum`
sample_data(ps)$mora=otu$`Moraxella sp.`
sample_data(ps)$m_lu=otu$`Micrococcus luteus`
sample_data(ps)$diet=otu$`Dietzia sp.`
sample_data(ps)$koc=otu$'Kocuria palustris'

sample_data(ps)$village=factor(sample_data(ps)$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))


data=data.frame(vill=sample_data(ps)$village, c_acne=sample_data(ps)$c_acnes,
s_epi=sample_data(ps)$s_epi, c_tuber=sample_data(ps)$c_tuber,
diet=sample_data(ps)$diet, koc=sample_data(ps)$koc,
mora=sample_data(ps)$mora,  mlu=sample_data(ps)$m_lu,site=sample_data(ps)$site_specific)

data2=filter(data, site == "Fh")

data2= data2 %>% select(1,2,8)


n=plot_correlation(
  data2,
  type = c("all"),
  maxcat = 20L,
  cor_args = list(),
  geom_text_args = list(),
  title = NULL,
  ggtheme = theme_gray(),
  theme_config = list(legend.position = "bottom", axis.text.x = element_text(angle = 90))
)

n=data.frame(n$data)
n2=dplyr::filter(n, Var1 %in% c("c_acne","mlu"))
n2=dplyr::filter(n2, !Var2 %in% c("c_acne","mlu"))


n2$Var2=gsub("vill_","",n2$Var2)

n2$Var2=factor(n2$Var2,levels=c("Kuala.Lumpur","Washington.DC","Village.F",
                                        "Village.D","Village.G","Village.AB",
                                        "Village.C"))
levels(n2$Var2)=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C")

n2$Var1=recode(n2$Var1, "c_acne"="C. acnes","m_lu"="S. epidermidis",
"c_tuber"="C. tuberculostearicum","mora"="Moraxella sp.",
"mlu"="M. luteus")

x<-ggplot(n2, aes(x=Var1,y = factor(Var2),fill=value))  +
  geom_tile(aes(fill= value)) +theme_classic()+  geom_text(aes(label = round(value, 2)), color = "black", size = 3, hjust = 0.5, vjust = 0.5)+
  labs(x="Species", y="Location", 
       fill="Pearson correlation coefficients")+
  theme(axis.text.y=element_text(size=8,face="bold"), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("royalblue", "lightblue1","white","red","darkred"),
                       breaks=c(-0.50,0,0.50),labels=c("-0.50","0","0.50"),
                       limits=c(-0.50,0.50))

a<- x+
  theme(axis.text.y = element_text(size=8,face="bold.italic"),
        axis.text.x =  element_text(size=8,face="bold",angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))

a=a+coord_flip()


#vf

data2=filter(data, site == "Vf")

data2= data2 %>% select(1,2,5)


n=plot_correlation(
  data2,
  type = c("all"),
  maxcat = 20L,
  cor_args = list(),
  geom_text_args = list(),
  title = NULL,
  ggtheme = theme_gray(),
  theme_config = list(legend.position = "bottom", axis.text.x = element_text(angle = 90))
)

n=data.frame(n$data)
n2=filter(n,!Var2 %in% c("c_acne","s_epi","c_tuber","mora","diet"))
n2=filter(n2, Var1 %in% c("c_acne","diet"))

n2$Var2=gsub("vill_","",n2$Var2)


n2$Var2=factor(n2$Var2,levels=c("Kuala.Lumpur","Washington.DC","Village.F",
                                        "Village.D","Village.G","Village.AB",
                                        "Village.C"))
levels(n2$Var2)=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C")

n2$Var1=recode(n2$Var1, "c_acne"="C. acnes","diet"="Dietzia sp.",
"c_tuber"="C. tuberculostearicum","mora"="Moraxella sp.",
"mlu"="M. luteus")

x<-ggplot(n2, aes(x=Var1,y = factor(Var2),fill=value)) + 
  geom_tile(aes(fill= value)) +theme_classic()+  geom_text(aes(label = round(value, 2)), color = "black", size = 3, hjust = 0.5, vjust = 0.5) +
  labs(x="Species", y="Location", 
       fill="Pearson correlation coefficients")+
  theme(axis.text.y=element_text(size=8,face="bold"), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("royalblue", "lightblue1","white", "red","darkred"),
                       breaks=c(-0.60,0,0.60),labels=c("-0.60","0","0.60"),
                       limits=c(-0.60,0.60))

b<- x+
  theme(axis.text.y = element_text(size=8,face="bold.italic"),
        axis.text.x =  element_text(size=8,face="bold",angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))

b=b+coord_flip()

#tw



data2=filter(data, site == "Tw")

data2= data2 %>% select(1,3,4,6,7)



n=plot_correlation(
  data2,
  type = c("all"),
  maxcat = 20L,
  cor_args = list(),
  geom_text_args = list(),
  title = NULL,
  ggtheme = theme_gray(),
  theme_config = list(legend.position = "bottom", axis.text.x = element_text(angle = 90))
)

n=data.frame(n$data)
n2=filter(n,!Var2 %in% c("c_acne","s_epi","c_tuber","mora","koc"))
n2=filter(n2, Var1 %in% c("s_epi","c_tuber","mora","koc"))

#n2=filter(n2, Var1 %in% c("c_acne","s_epi","c_tuber","mora","mlu"))

n2$Var2=gsub("vill_","",n2$Var2)

n2$Var2=factor(n2$Var2,levels=c("Kuala.Lumpur","Washington.DC","Village.F",
                                        "Village.D","Village.G","Village.AB",
                                        "Village.C"))
levels(n2$Var2)=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C")


n2$Var1=recode(n2$Var1, "c_acne"="C. acnes","s_epi"="S. epidermidis",
"c_tuber"="C. tuberculostearicum","mora"="Moraxella sp.",
"koc"="K. palustris")

x<-ggplot(n2, aes(x=Var1,y = factor(Var2),fill=value)) + 
  geom_tile(aes(fill= value)) +theme_classic() +  geom_text(aes(label = round(value, 2)), color = "black", size = 3, hjust = 0.5, vjust = 0.5)+
  labs(x="Species", y="Location", 
       fill="Pearson correlation coefficients")+
  theme(axis.text.y=element_text(size=8,face="bold"), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("royalblue", "lightblue1","white", "red","darkred"),
                       breaks=c(-0.56,0,0.56),labels=c("-0.56","0","0.56"),
                       limits=c(-0.56,0.56))

c<- x+
  theme(axis.text.y = element_text(size=8,face="bold.italic"),
        axis.text.x =  element_text(size=8,face="bold",angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))

c=c+coord_flip()


q=ggarrange(a+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none"),
b+theme(legend.position="none"),heights=c(1,1.6),nrow=2,ncol=1)

q2=ggarrange(q,c+theme(legend.position="none"),nrow=1,ncol=2)

ggsave("Fig_S4A-B.svg",q2,height=2.6,width=7.7)
