library(pacman)

pacman::p_load(readr,plyr, dplyr, ggpubr, vegan, plotly, plyr, ggplot2)


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")



ps=readRDS("phyloseq_0505.rda")


ps=subset_samples(ps,!village %in%  c("Washington DC","Kuala Lumpur"))
ps=subset_samples(ps,site_specific != "Ctrl")
ps=subset_samples(ps,Disease %in% c("Negative","Tinea imbricata"))

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

 


n=data.frame(library_id=sample_data(ps)$library_id,au=sample_data(ps)$a_u,
subject_id=sample_data(ps)$subject_id, Disease=sample_data(ps)$Disease,site=sample_data(ps)$site_specific)


otu2$library_id=rownames(otu2)
df=left_join(otu2,n)

data4=melt(df)


data4$Disease=factor(data4$Disease,
                     levels=c("Negative","Tinea imbricata"))
data4$site=factor(data4$site,levels=c("Fh","Vf","Tw"))
levels(data4$site)=c("Forehead","Volar forearm","Toeweb")


data4$group_status <- ifelse(data4$Disease == "Negative", "Healthy control", ifelse(data4$Disease == "Tinea imbricata", data4$au,NA))


data4$group_status=recode(data4$group_status,"A"="Lesional skin","U"="Nonlesional skin")
library(pacman)

data4$group_status=factor(data4$group_status,levels=c("Healthy control","Nonlesional skin","Lesional skin"))

data4$value_l=log(1+data4$value)

data4=filter(data4, variable %in% c("Staphylococcus.caprae",
"Staphylococcus.warneri","Corynebacterium.flavescens","Corynebacterium.terpenotabidum"))


data4$value_l=log(1+data4$value)

data4=filter(data4, site == "Forehead")

data4$variable=factor(data4$variable,levels=c("Staphylococcus.caprae",
"Staphylococcus.warneri","Corynebacterium.flavescens","Corynebacterium.terpenotabidum"))


comparison=list(c("Healthy control","Nonlesional skin"),c("Nonlesional skin","Lesional skin"),
                c("Healthy control","Lesional skin"))



a <- ggboxplot(data4, x ="group_status" , y = "value_l",
               color = "group_status"               )+theme_classic()+
  labs(x="Status",y="Abundance (log normalized)",
       color="Status")+
scale_color_manual(values=c("royalblue","salmon","red3"))+
facet_wrap(~ variable,nrow=1) + stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Healthy control",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x =element_blank(),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size=8, face="italic"),
         plot.title=element_text(size=8),
strip.background = element_rect(fill="white"))




data4=melt(df)


data4$Disease=factor(data4$Disease,
                     levels=c("Negative","Tinea imbricata"))
data4$site=factor(data4$site,levels=c("Fh","Vf","Tw"))
levels(data4$site)=c("Forehead","Volar forearm","Toeweb")


data4$group_status <- ifelse(data4$Disease == "Negative", "Healthy control", ifelse(data4$Disease == "Tinea imbricata", data4$au,NA))


data4$group_status=recode(data4$group_status,"A"="Lesional skin","U"="Nonlesional skin")
library(pacman)

data4$group_status=factor(data4$group_status,levels=c("Healthy control","Nonlesional skin","Lesional skin"))

data4$value_l=log(1+data4$value)

data4=filter(data4, variable %in% c("Staphylococcus.caprae",
"Brucella.anthropi","Corynebacterium.flavescens","Corynebacterium.terpenotabidum"))


data4$value_l=log(1+data4$value)

data4=filter(data4, site == "Volar forearm")

data4$variable=factor(data4$variable,levels=c("Staphylococcus.caprae",
"Brucella.anthropi","Corynebacterium.flavescens","Corynebacterium.terpenotabidum"))


comparison=list(c("Healthy control","Nonlesional skin"),c("Nonlesional skin","Lesional skin"),
                c("Healthy control","Lesional skin"))



b <- ggboxplot(data4, x ="group_status" , y = "value_l",
               color = "group_status"              )+theme_classic()+
  labs(x="Status",y="Abundance (log normalized)",
       color="Status")+
scale_color_manual(values=c("royalblue","salmon","red3"))+
facet_wrap(~ variable,nrow=1) + stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Healthy control",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x =element_blank(),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size=8, face="italic"),
         plot.title=element_text(size=8),
strip.background = element_rect(fill="white"))



data4=melt(df)


data4$Disease=factor(data4$Disease,
                     levels=c("Negative","Tinea imbricata"))
data4$site=factor(data4$site,levels=c("Fh","Vf","Tw"))
levels(data4$site)=c("Forehead","Volar forearm","Toeweb")


data4$group_status <- ifelse(data4$Disease == "Negative", "Healthy control", ifelse(data4$Disease == "Tinea imbricata", data4$au,NA))


data4$group_status=recode(data4$group_status,"A"="Lesional skin","U"="Nonlesional skin")
library(pacman)

data4$group_status=factor(data4$group_status,levels=c("Healthy control","Nonlesional skin","Lesional skin"))

data4$value_l=log(1+data4$value)

data4=filter(data4, variable %in% c("Kytococcus.sp."))

data4$value_l=log(1+data4$value)

data4=filter(data4, site == "Toeweb")

data4$variable=factor(data4$variable,levels=c("Kytococcus.sp."))


comparison=list(c("Healthy control","Nonlesional skin"),c("Nonlesional skin","Lesional skin"),
                c("Healthy control","Lesional skin"))



c <- ggboxplot(data4, x ="group_status" , y = "value_l",
               color = "group_status"              )+theme_classic()+
  labs(x="Status",y="Abundance (log normalized)",
       color="Status")+
scale_color_manual(values=c("royalblue","salmon","red3"))+
facet_wrap(~ variable,nrow=1) + stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "Healthy control",angle=90)+
  theme( axis.text.x = element_blank(),
         axis.text.y = element_text(size=8),
         axis.title.x =element_blank(),
         axis.title.y = element_text(size=8),
         strip.text.x = element_text(size=8, face="italic"),
         plot.title=element_text(size=8),
strip.background = element_rect(fill="white"))




q=ggarrange(a+ylim(0,10)+theme(legend.position="none"),b+
ylim(0,10)+theme(legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank()),
c+theme(legend.position="none",axis.text.y=element_blank(),
axis.title.y=element_blank()),nrow=1,ncol=3, widths=c(1,0.85,0.2))

ggsave("Fig_S6E.svg",q,dpi=320,height=2,width=7.7)

