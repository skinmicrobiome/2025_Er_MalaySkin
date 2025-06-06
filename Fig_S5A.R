library(pacman)

pacman::p_load(tidyverse, dplyr, ggplot2,ggsci,vegan,reshape2,phyloseq,ggpubr,readr,tidyverse)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db/input_file_r_code_coded")

ps=readRDS("phyloseq_S5A.rda")

ps = tax_glom(ps, "Species")




ps1=subset_samples(ps, village =="Washington DC")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")


components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")





a= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(49.26%)",
y="PC2 (22.07%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))
#lumpur 


ps1=subset_samples(ps, village =="Kuala Lumpur")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)
summary(prin_comp)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")


components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")





b= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(60.83%)",
y="PC2 (13.31%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))
#leg


ps1=subset_samples(ps, village =="Village F")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)
summary(prin_comp)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")



components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toe web")


summary(prin_comp)


c= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(61.29%)",
y="PC2 (16.58%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))


#garam

ps1=subset_samples(ps, village =="Village D")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)
summary(prin_comp)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")



components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toe web")

summary(prin_comp)




d= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(47.13%)",
y="PC2 (15.08%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))
#temen


ps1=subset_samples(ps, village =="Village G")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)
summary(prin_comp)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")



components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")

summary(prin_comp)



e= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(71.06%)",
y="PC2 (7.43%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))


#kuala koh



ps1=subset_samples(ps, village =="Village AB")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)
summary(prin_comp)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")



components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toe web")


summary(prin_comp)


f= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(51.42%)",
y="PC2 (19.32%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))

#tn




ps1=subset_samples(ps, village =="Village C")

otu=as.data.frame(t(otu_table(ps1)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)
summary(prin_comp)
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps1)$village,sample_data(ps1)$site_specific,
                   sample_data(ps1)$Disease)
colnames(components)=c("PC1","PC2","PC3","village","site","disease")

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")



summary(prin_comp)

g= ggplot(components, aes(x = PC1, y = PC2, color = site)) +
  geom_point( alpha=0.70,size=1.5) +labs(color="Body site",x="PC1(53.90%)",
y="PC2 (20.01%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("seagreen","orange3","purple3"))  +labs(color="Body site")+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  guides(shape = guide_legend(override.aes = list(size=8)))

q=ggarrange(a+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),
b+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),
c+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),
d+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),
e+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),
f+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),
g+theme_classic()+theme(axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
legend.position="none"),nrow=2,ncol=4,align="hv")


ggsave("Fig_S5A.svg",q,height=3,width=8.5)

ggsave("Fig_S5_legend.svg",get_legend(c),height=2,width=2)



