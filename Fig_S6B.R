library(pacman)

pacman::p_load(tidyverse, dplyr, ggplot2,ggsci,vegan,reshape2,phyloseq,ggpubr,readr,tidyverse)


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


ps2=readRDS("phyloseq_S6B.rda")

ps3=subset_samples(ps2, site_specific %in% c("Fh"))


#ps2=subset_samples(ps1, site_specific =="Fh")
metadata<-sample_data(ps3)

otu <- as.data.frame(t(otu_table(ps3)))

dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]

summary(prin_comp)

components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps3)$village,sample_data(ps3)$site_specific,
                   sample_data(ps3)$Disease, sample_data(ps3)$a_u)
colnames(components)=c("PC1","PC2","PC3","village","site","disease","a_u")



components$village=factor(components$village,levels=c("Bethesda","Kuala Lumpur","Lubuk Legong",
                                                      "Garam","Tementong","Kuala Koh","Taman Negara"))

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")

components$disease=recode(components$disease,"Negative"="Healthy control")

components$disease=factor(components$disease,levels=c("Healthy control", "Tinea imbricata"
                                                                                                       ))

components$a_u=factor(components$a_u,levels=c("A","U"))
levels(components$a_u)=c("With lesion","Without lesion")



components$group_status <- ifelse(components$disease == "Healthy control", "Healthy control", ifelse(components$disease == "Tinea imbricata", components$a_u,NA))


components$group_status=recode(components$group_status,"1"="Lesional  skin","2"="Nonlesional skin")


components$group_status=factor(components$group_status, levels=c("Healthy control","Nonlesional skin","Lesional  skin"))


a= ggplot(components, aes(x = PC1, y = PC2, color = group_status)) +
  geom_point( alpha=0.70,size=3) +labs(color="Status",x="PC1(55.96%)",
y="PC2 (11.94%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("royalblue","salmon2","violetred3")) +labs(color="Status")+
##facet_grid(~site)+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=5)))

#vf

ps3=subset_samples(ps2, site_specific %in% c("Vf"))


metadata<-sample_data(ps3)

otu <- as.data.frame(t(otu_table(ps3)))

dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]

summary(prin_comp)

components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps3)$village,sample_data(ps3)$site_specific,
                   sample_data(ps3)$Disease, sample_data(ps3)$a_u)
colnames(components)=c("PC1","PC2","PC3","village","site","disease","a_u")



components$village=factor(components$village,levels=c("Bethesda","Kuala Lumpur","Lubuk Legong",
                                                      "Garam","Tementong","Kuala Koh","Taman Negara"))

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")

components$disease=recode(components$disease,"Negative"="Healthy control")

components$disease=factor(components$disease,levels=c("Healthy control", "Tinea imbricata"
                                                                                                       ))

components$a_u=factor(components$a_u,levels=c("A","U"))
levels(components$a_u)=c("With lesion","Without lesion")



components$group_status <- ifelse(components$disease == "Healthy control", "Healthy control", ifelse(components$disease == "Tinea imbricata", components$a_u,NA))


components$group_status=recode(components$group_status,"1"="Lesional  skin","2"="Nonlesional skin")


components$group_status=factor(components$group_status, levels=c("Healthy control","Nonlesional skin","Lesional  skin"))


b= ggplot(components, aes(x = PC1, y = PC2, color = group_status)) +
  geom_point( alpha=0.70,size=3) +labs(color="Status",x="PC1(53.62%)",
y="PC2 (10.16%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("royalblue","salmon2","violetred3")) +labs(color="Status")+
##facet_grid(~site)+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=5)))

#tw


ps3=subset_samples(ps2, site_specific %in% c("Tw"))


metadata<-sample_data(ps3)

otu <- as.data.frame(t(otu_table(ps3)))

dist <- vegdist(decostand(otu, "hellinger"), "bray")


prin_comp <- prcomp(dist, rank. = 3)

components <- prin_comp[["x"]]

summary(prin_comp)

components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, sample_data(ps3)$village,sample_data(ps3)$site_specific,
                   sample_data(ps3)$Disease, sample_data(ps3)$a_u)
colnames(components)=c("PC1","PC2","PC3","village","site","disease","a_u")



components$village=factor(components$village,levels=c("Bethesda","Kuala Lumpur","Lubuk Legong",
                                                      "Garam","Tementong","Kuala Koh","Taman Negara"))

components$site=factor(components$site,levels=c("Fh","Vf","Tw"))
levels(components$site)=c("Forehead","Volar forearm","Toeweb")

components$disease=recode(components$disease,"Negative"="Healthy control")

components$disease=factor(components$disease,levels=c("Healthy control", "Tinea imbricata"
                                                                                                       ))

components$a_u=factor(components$a_u,levels=c("A","U"))
levels(components$a_u)=c("With lesion","Without lesion")


components$group_status <- ifelse(components$disease == "Healthy control", "Healthy control", ifelse(components$disease == "Tinea imbricata", components$a_u,NA))


components$group_status=recode(components$group_status,"1"="Lesional  skin","2"="Nonlesional skin")


components$group_status=factor(components$group_status, levels=c("Healthy control","Nonlesional skin","Lesional  skin"))


c= ggplot(components, aes(x = PC1, y = PC2, color = group_status)) +
  geom_point( alpha=0.70,size=3) +labs(color="Status",x="PC1(67.90%)",
y="PC2 (9.62%)")+
  theme_classic() + theme(legend.position="bottom")+
  scale_colour_manual(values=c("royalblue","salmon2","violetred3")) +labs(color="Status")+
##facet_grid(~site)+
  theme( axis.text.x = element_text(size=8,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         strip.text.x = element_text(size =8, face="bold"),
         plot.title = element_text(size=8,face="bold"),
         legend.text = element_text(size=8,face="bold"),
         legend.title = element_text(size=8,face="bold"),
         legend.position = "bottom")+
  guides(colour = guide_legend(override.aes = list(size=5)))

q=ggarrange(a+theme(legend.position="none"),b+theme(legend.position="none"),
c+theme(legend.position="none"),nrow=1,ncol=3)


q3=ggarrange(q2+theme(legend.position="none"),NULL, nrow=2,ncol=1,
heights=c(1,0.17))

q3=ggarrange(q,q3, widths=c(1,0.4))

ggsave(file="merge_1_0605.svg", plot=q3, width=7.7, height=2,dpi=320)



q2=ggarrange(q,get_legend(c),nrow=2,ncol=1,heights=c(1,0.1))

ggsave(file="beta_ti_det_pcoa_0407.svg", plot=q2, width=5.5, height=3,dpi=320)


ggsave(file="beta_ti_det_pcoa_2002_leg.svg", plot=get_legend(c), width=5, height=2,dpi=320)

