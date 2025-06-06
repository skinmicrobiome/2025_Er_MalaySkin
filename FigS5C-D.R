library(pacman)

pacman::p_load(readr,plyr, dplyr, ggpubr, vegan, plotly, plyr, ggplot2)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db	")

z2=readr::read_tsv("vill_150324/significant_results.tsv")

z2$Sig_assoc=-log(z2$qval)*sign(z2$coef)

z2$abs_sig=abs(z2$Sig_assoc)


z3=filter(z2,value =="Village C")
z3=z3[order(-z3$abs_sig),]
z3=head(z3,n=10)


z4=filter(z2,value =="Village AB")
z4=z4[order(-z4$abs_sig),]
z4=head(z4,n=10)

z7=filter(z2,value =="Village F")
z7=z7[order(-z7$abs_sig),]
z7=head(z7,n=10)


z8=filter(z2,value =="Washington DC")
z8=z8[order(-z8$abs_sig),]
z8=head(z8,n=10)


library(dplyr)

z3 <- z3 %>%
  group_by(value) %>%
  arrange(desc(abs_sig)) %>%
  ungroup()

ps=readRDS("phyloseq_0505.rda")

ps=subset_samples(ps,site_specific != "Ctrl")
ps=microbiome::transform(ps,"rank")
tax=data.frame(tax_table(ps))
tax$sp2=paste(tax$Genus,tax$Species)
tax2=tax_table(tax)
taxa_names(tax2)=rownames(tax)
tax_table(ps)=tax2
colnames(tax_table(ps))=c("Kingdon","Phylum","Class","Order","Family",
                          "Genus","Species","sp2")

otu=data.frame(otu_table(ps))

rownames(otu)=tax$sp2
otu2=data.frame(t(otu))

 


n=data.frame(library_id=sample_data(ps)$library_id,
subject_id=sample_data(ps)$subject_id, village=sample_data(ps)$village,site=sample_data(ps)$site_specific)


otu2$library_id=rownames(otu2)
df=dplyr::left_join(otu2,n)

df3=melt(df)

df4=filter(df3, variable %in% z3$feature)



df4$site=factor(df4$site,levels=c("Fh","Vf","Tw"))
levels(df4$site)=c("Forehead","Volar forearm","Toeweb")


df5=filter(df4, village %in% c("Kuala Lumpur","Village C"))

x<-ggplot(df5, aes(x=village,y = factor(variable, levels=unique(z3$feature)),fill=value)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Village", y="Species", 
       fill="Significant associations")+
  theme(axis.text.y=element_text(size=8,face="italic"), 
        legend.position = "bottom")+facet_grid(~site)+
  scale_fill_gradientn(colours = c("white", "orange3","violetred3"),
                       breaks=c(0,6000),labels=c("0",
                                                            "6000"),
                       limits=c(0,6000))

a<- x+
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=0.5,size=8),
        axis.title.x = element_text(size=8),
        strip.text.x = element_text(size=7.5),
        axis.title.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "right",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))

#kk

df4=filter(df3, variable %in% z4$feature)


df4$site=factor(df4$site,levels=c("Fh","Vf","Tw"))
levels(df4$site)=c("Forehead","Volar forearm","Toeweb")


df5=filter(df4, village %in% c("Kuala Lumpur","Village AB"))

x<-ggplot(df5, aes(x=village,y = factor(variable, levels=unique(z4$feature)),fill=value)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Village", y="Species", 
       fill="Significant associations")+
  theme(axis.text.y=element_text(size=8,face="italic"), 
        legend.position = "bottom")+facet_grid(~site)+
  scale_fill_gradientn(colours = c("white", "orange3","violetred3"),
                       breaks=c(0,6000),labels=c("0",
                                                            "6000"),
                       limits=c(0,6000))
b<- x+
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=0.5,size=8),
        axis.title.x = element_text(size=8),
        strip.text.x = element_text(size=7.5),
        axis.title.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "right",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))
#ll


df4=filter(df3, variable %in% z7$feature)



df4$site=factor(df4$site,levels=c("Fh","Vf","Tw"))
levels(df4$site)=c("Forehead","Volar forearm","Toeweb")


df5=filter(df4, village %in% c("Kuala Lumpur","Village F"))

x<-ggplot(df5, aes(x=village,y = factor(variable),fill=value)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Village", y="Species", 
       fill="Significant associations")+
  theme(axis.text.y=element_text(size=8,face="italic"), 
        legend.position = "bottom")+facet_grid(~site)+
  scale_fill_gradientn(colours = c("white", "orange3","violetred3"),
                       breaks=c(0,6000),labels=c("0",
                                                            "6000"),
                       limits=c(0,6000))
e<- x+
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=0.5,size=8),
        axis.title.x = element_text(size=8),
        strip.text.x = element_text(size=7.5),
        axis.title.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "right",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))

#be

df4=filter(df3, variable %in% z8$feature)


df4$site=factor(df4$site,levels=c("Fh","Vf","Tw"))
levels(df4$site)=c("Forehead","Volar forearm","Toeweb")


df5=filter(df4, village %in% c("Kuala Lumpur","Washington DC"))

x<-ggplot(df5, aes(x=village,y = factor(variable, levels=unique(z8$feature)),fill=value)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Village", y="Species", 
       fill="Significant associations")+
  theme(axis.text.y=element_text(size=8,face="italic"), 
        legend.position = "bottom")+facet_grid(~site)+
  scale_fill_gradientn(colours = c("white", "orange3","violetred3"),
                       breaks=c(0,6000),labels=c("0",
                                                            "6000"),
                       limits=c(0,6000))
f<- x+
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=0.5,size=8),
        axis.title.x = element_text(size=8),
        strip.text.x = element_text(size=7.5),
        axis.title.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "right",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))



q=ggarrange(a+theme(legend.position="none"),b+theme(legend.position="none"),e+theme(legend.position="none"),
f+theme(legend.position="none"),nrow=2,ncol=2)

a=a+theme(legend.position="bottom")


q2=ggarrange(q,get_legend(a),nrow=2,ncol=1, heights=c(1,0.1))

ggsave("Fig_S5C-D.svg",q2,height=4,width=8)





