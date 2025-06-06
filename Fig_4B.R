library(pacman)

pacman::p_load(tidyverse, dplyr,DT, mia, scater,  phyloseq, processx)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

ps=readRDS("rare_1000_spkm_1803.rda")


ps=readRDS("phyloseq_0505.rda")

ps = tax_glom(ps, "Species")

ps1=subset_samples(ps,site_specific == "Fh")

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



components$village=factor(components$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")


colors = c("#FFBC00","#FF8000","#b1ffb1",
                          "#27ff27","#00b100","#25A18E","#004e00")


summary(prin_comp)


a= ggplot(components, aes(x = PC1, y = PC2, color = village)) +
  geom_point( alpha=0.70,size=3) +labs(color="Location",x="PC1(65.99%)",
y="PC2 (14.69%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=colors)  +labs(color="Location",title="Forehead")+
  theme( axis.text.x = element_text(size=8),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8),
         legend.text = element_text(size=8),
         legend.title = element_text(size=8),
         legend.position = "right")



#vf

ps1=subset_samples(ps,site_specific == "Vf")

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



components$village=factor(components$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                                      "Village D","Village G","Village AB","Village C"))

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")

summary(prin_comp)

b= ggplot(components, aes(x = PC1, y = PC2, color = village)) +
  geom_point( alpha=0.70,size=3) +labs(color="Location",x="PC1(65.45%)",
y="PC2 (8.80%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=colors)  +labs(color="Location",title="Volar forearm")+
  theme( axis.text.x = element_text(size=8),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8),
         legend.text = element_text(size=8),
         legend.title = element_text(size=8),
         legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size=5)))


#tw


ps1=subset_samples(ps,site_specific == "Tw")

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



components$village=factor(components$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                                      "Village D","Village G","Village AB","Village C"))


components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")

summary(prin_comp)

c= ggplot(components, aes(x = PC1, y = PC2, color = village)) +
  geom_point( alpha=0.70,size=3) +labs(color="Location",x="PC1(72.39%)",
y="PC2 (9.70%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=colors)  +labs(color="Location",title="Toeweb")+
  theme( axis.text.x = element_text(size=8),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8),
         legend.text = element_text(size=8),
         legend.title = element_text(size=8),
         legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size=5)))


q2=ggarrange(a+theme(legend.position="none"),
b+theme(legend.position="none"),
c+theme(legend.position="none"),nrow=1,ncol=3)

c=c+theme(legend.position="bottom")+
 guides(color = guide_legend(override.aes = list(size=5),nrow=1))



q2=ggarrange(q2,get_legend(c),heights=c(1,0.2),nrow=2,ncol=1)


ggsave(file="Fig_4B.svg", plot=q2, width=6, height=1.5,dpi=320)



