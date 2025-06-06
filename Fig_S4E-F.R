
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
colnames(tax_table(ps))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species","sp2")



ps2=subset_taxa(ps, sp2 %in% c("Moraxella sp.","Micrococcus luteus","Corynebacterium tuberculostearicum",
"Cutibacterium acnes","Staphylococcus epidermidis","Dietzia sp.","Kocuria palustris"
))


ps2=microbiome::transform(ps2,"Z")

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




ps1 = tax_glom(ps, "sp2")

ps1


ps1=subset_samples(ps1,site_specific != "Ctrl")

ps2=subset_samples(ps1,site_specific == "Fh")


metadata<-sample_data(ps2)

otu <- as.data.frame(t(otu_table(ps2)))

dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]


summary(prin_comp)


components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps2)$village,sample_data(ps2)$site_specific,
                   sample_data(ps2)$Disease,sample_data(ps2)$c_acnes,
sample_data(ps2)$s_epi,sample_data(ps2)$c_tuber,sample_data(ps2)$mora,sample_data(ps2)$m_lu,
sample_data(ps2)$diet,sample_data(ps2)$koc)
colnames(components)=c("PC1","PC2","PC3","village","site","disease","c_acnes",
"s_epi","c_tuber","mora","m_lu","diet","koc")


tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")


components$OAUC=recode(components$village,"Washington DC"="Urban","Kuala Lumpur"="Urban",.default="Orang Asli")
components$OAUC=factor(components$OAUC,levels=c("Urban","Orang Asli"))

a= ggplot(components, aes(x = PC1, y = PC2, color=c_acnes,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="C. acnes abundance (log normalized",
shape="Urbanization scores",x="PC1(64.82%)",
y="PC2 (14.59%)",title="Forehead")+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-2.5,2.5),labels=c("0","18"
                        ),
                        limits=c(-2.5,2.5))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")
  


b= ggplot(components, aes(x = PC1, y = PC2, color=m_lu,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="M. luteus abundance (log normalized",
shape="Urbanization scores",x="PC1(64.82%)",
y="PC2 (14.59%)",title="Forehead")+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-3.5,3.5),labels=c("-3.5","3.5"
                        ),
                        limits=c(-3.5,3.5))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")
  






#vf



ps2=subset_samples(ps1,site_specific == "Vf")


metadata<-sample_data(ps2)

otu <- as.data.frame(t(otu_table(ps2)))

dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]


summary(prin_comp)


components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps2)$village,sample_data(ps2)$site_specific,
                   sample_data(ps2)$Disease,sample_data(ps2)$c_acnes,
sample_data(ps2)$s_epi,sample_data(ps2)$c_tuber,sample_data(ps2)$mora,sample_data(ps2)$m_lu,
sample_data(ps2)$diet,sample_data(ps2)$koc)
colnames(components)=c("PC1","PC2","PC3","village","site","disease","c_acnes",
"s_epi","c_tuber","mora","m_lu","diet","koc")



tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")


components$OAUC=recode(components$village,"Washington DC"="Urban","Kuala Lumpur"="Urban",.default="Orang Asli")
components$OAUC=factor(components$OAUC, levels=c("Urban","Orang Asli"))

c= ggplot(components, aes(x = PC1, y = PC2, color=c_acnes,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="C. acnes abundance (log normalized",
shape="Urbanization scores",x="PC1(64.88%)",
y="PC2 (8.83%)",title="Volar forearm")+ scale_shape_manual(values=c(16, 17))+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-1.7,1.7),labels=c("-1.7","1.7"
                        ),
                        limits=c(-1.7,1.7))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")

d= ggplot(components, aes(x = PC1, y = PC2, color=diet,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="Dietzia sp. abundance (log normalized",
shape="Urbanization scores",x="PC1(64.88%)",
y="PC2 (8.83%)",title="Volar forearm")+ scale_shape_manual(values=c(16, 17))+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-1.8,1.8),labels=c("-1.8","1.8"
                        ),
                        limits=c(-1.8,1.8))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")
  

#tw



ps2=subset_samples(ps1,site_specific == "Tw")


metadata<-sample_data(ps2)

otu <- as.data.frame(t(otu_table(ps2)))

dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]


summary(prin_comp)


components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps2)$village,sample_data(ps2)$site_specific,
                   sample_data(ps2)$Disease,sample_data(ps2)$c_acnes,
sample_data(ps2)$s_epi,sample_data(ps2)$c_tuber,sample_data(ps2)$mora,sample_data(ps2)$m_lu,
sample_data(ps2)$diet,sample_data(ps2)$koc)
colnames(components)=c("PC1","PC2","PC3","village","site","disease","c_acnes",
"s_epi","c_tuber","mora","m_lu","diet","koc")



tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")



components$OAUC=recode(components$village,"Washington DC"="Urban","Kuala Lumpur"="Urban",.default="Orang Asli")
components$OAUC=factor(components$OAUC, levels=c("Orang Asli","Urban"))


e= ggplot(components, aes(x = PC1, y = PC2, color=s_epi,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="S. epidermidis abundance (log normalized",
shape="Urbanization scores",x="PC1(72.50%)",
y="PC2 (9.55%)",title="Toeweb")+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-3.1,3.1),labels=c("-3.1","3.1"
                        ),
                        limits=c(-3.1,3.1))+ scale_shape_manual(values=c(17, 16))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")

f= ggplot(components, aes(x = PC1, y = PC2, color=c_tuber,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="C. tuberculostearicum abundance (log normalized",
shape="Urbanization scores",x="PC1(72.50%)",
y="PC2 (9.55%)",title="Toeweb")+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-3,3),labels=c("-3","3"
                        ),
                        limits=c(-3,3))+ scale_shape_manual(values=c(17, 16))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")


g= ggplot(components, aes(x = PC1, y = PC2, color=mora,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="Moraxella sp. abundance (log normalized",
shape="Urbanization scores",x="PC1(72.50%)",
y="PC2 (9.55%)",title="Toeweb")+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-2.7,2.7),labels=c("-2.7","2.7"
                        ),
                        limits=c(-2.7,2.7))+ scale_shape_manual(values=c(17, 16))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = "none")

h= ggplot(components, aes(x = PC1, y = PC2, color=koc,shape = as.factor(OAUC))) +
  geom_point( alpha=0.70,size=3) +labs(color="K. palustris abundance (log normalized",
shape="Urbanization scores",x="PC1(72.50%)",
y="PC2 (9.55%)",title="Toeweb")+
  theme_classic() + theme(legend.position="bottom") +
  scale_color_gradientn(colours = c("white","pink","darkred"),
                        breaks=c(-2.6,2.6),labels=c("-2.6","2.6"
                        ),
                        limits=c(-2.6,2.6))+ scale_shape_manual(values=c(17, 16))+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(shape = guide_legend(override.aes = list(size=5)))
  



q=ggarrange(a+theme(legend.position="none"),c+theme(legend.position="none"),
e+theme(legend.position="none"),f+theme(legend.position="none"),nrow=1,ncol=4)


ggsave("Fig_S4E_0407.svg",q,height=2,width=8)


a=a+theme(legend.position="bottom")
c=c+theme(legend.position="bottom")
e=e+theme(legend.position="bottom")
f=f+theme(legend.position="bottom")

q2=ggarrange(get_legend(a),get_legend(c),get_legend(e),get_legend(f),
nrow=4,ncol=1)


ggsave("Fig_S4E_legend_1703.svg",q2,height=7,width=7)





q3=ggarrange(b+theme(legend.position="none"),d+theme(legend.position="none"),
g+theme(legend.position="none"),h+theme(legend.position="none"),nrow=1,ncol=4)

ggsave("Fig_S4F_0407.svg",q3,height=2,width=8)



b=b+theme(legend.position="bottom")
d=d+theme(legend.position="bottom")
g=g+theme(legend.position="bottom")
h=h+theme(legend.position="bottom")

q4=ggarrange(get_legend(b),get_legend(d),get_legend(g),get_legend(h),
nrow=4,ncol=1)


ggsave("Fig_S4F_legend.svg",q4,height=7,width=7)



