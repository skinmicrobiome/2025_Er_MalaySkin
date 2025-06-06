library(pacman)

pacman::p_load(tidyverse, dplyr,DT, mia, scater,  phyloseq, vegan,microbiome)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")



ps=readRDS("phyloseq_0505.rda")


#OA_Only

ps=subset_samples(ps,!village %in%  c("Washington DC","Kuala Lumpur"))



ps

ps=subset_samples(ps, site_specific != "Ctrl")

meta=data.frame(sample_data(ps))


ps1 = tax_glom(ps, "Species")

Body_site=as.factor(sample_data(ps)$site_specific)
Village=as.factor(sample_data(ps)$village)
Disease=as.factor(sample_data(ps)$Disease)

BMI=as.factor(sample_data(ps)$Category)
Subtribe=as.factor(sample_data(ps)$Subtribe)
Shoes_type=as.factor(sample_data(ps)$shoes)
Age=as.factor(sample_data(ps)$Group)
Pets=as.factor(sample_data(ps)$Pets_1)
Smoking=as.factor(sample_data(ps)$Tobacco)
Gender=as.factor(sample_data(ps)$Gender)
Antifungals=as.factor(sample_data(ps)$Antifungal_3_months)
Presence_of_lesion=as.factor(sample_data(ps)$a_u)
Family_history=as.factor(sample_data(ps)$Family_skin_disease)
Drugs=as.factor(sample_data(ps)$Drug_3_months)
oilments=as.factor(sample_data(ps)$Oilment_cream_3_months)
Group=as.factor(sample_data(ps)$Group)
depth=as.numeric(log(sample_data(ps)$log_s))

otu=as.data.frame(t(otu_table(ps)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")

a=adonis2(dist~Body_site+Village+Disease+BMI+Subtribe+Shoes_type+Age+
Pets+Smoking+Gender+Antifungals+Family_history+Drugs+oilments, permutations=10000, by="margin")

a=data.frame(a)

colnames(a)=c("DF","SumOfSqs","R2","F","p_val")
a$q_val=p.adjust(a$p_val,method="BH")
write.csv(a,"adonis_oa_1703.csv")

#mdmr
library(MDMR)

z=otu_table(ps)

tax=data.frame(tax_table(ps))

rownames(z)=tax$Species
z=t(z)
z=data.matrix(z)

D <- dist(z, method = "euclidean")

X=data.frame(sample_data(ps1))
x=data.frame(Body_site,Village,Disease,BMI,Subtribe,Shoes_type,Age,
Pets,Smoking,Gender,Antifungals,Presence_of_lesion,Family_history,Drugs,oilments)

G <- gower(D)
lambda <- eigen(G, only.values = T)$values
mdmr.res2 <- mdmr(X =x, G = G, lambda = lambda)
b=summary(mdmr.res2)

b=data.frame(b)
colnames(b)=c("statistic","DF","R2","p_val")
b$q_val=p.adjust(b$p_val,method="BH")
write.csv(b,"mdmr_oa_1703.csv")


#all

library(pacman)

pacman::p_load(tidyverse, dplyr,DT, mia, scater,  phyloseq, processx)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

ps=readRDS("rare_1000_spkm_1803.rda")



ps

ps=subset_samples(ps, site_specific != "Ctrl")


meta=data.frame(sample_data(ps))



Body_site=as.factor(sample_data(ps)$site_specific)
Village=as.factor(sample_data(ps)$village)
Disease=as.factor(sample_data(ps)$Disease)

BMI=as.factor(sample_data(ps)$Category)
Subtribe=as.factor(sample_data(ps)$Subtribe)
Shoes_type=as.factor(sample_data(ps)$shoes)
Age=as.factor(sample_data(ps)$Group)
Pets=as.factor(sample_data(ps)$Pets_1)
Smoking=as.factor(sample_data(ps)$Tobacco)
Gender=as.factor(sample_data(ps)$Gender)
Antifungals=as.factor(sample_data(ps)$Antifungal_3_months)
Presence_of_lesion=as.factor(sample_data(ps)$a_u)
Family_history=as.factor(sample_data(ps)$Family_skin_disease)
Drugs=as.factor(sample_data(ps)$Drug_3_months)
oilments=as.factor(sample_data(ps)$Oilment_cream_3_months)
Group=as.factor(sample_data(ps)$Group)

otu=as.data.frame(t(otu_table(ps)))
dist <- vegdist(decostand(otu, "hellinger"), "bray")

a=adonis2(dist~Body_site+Village+Disease+BMI+Subtribe+Shoes_type+Age+
Pets+Smoking+Gender+Antifungals+Family_history+Drugs+oilments, permutations=10000, by="margin")

a=data.frame(a)

colnames(a)=c("DF","SumOfSqs","R2","F","p_val")
a$q_val=p.adjust(a$p_val,method="BH")
write.csv(a,"adonis_all_1703_1000.csv")

#mdmr
library(MDMR)

z=otu_table(ps)

tax=data.frame(tax_table(ps))

rownames(z)=tax$Species
z=t(z)
z=data.matrix(z)

D <- dist(z, method = "euclidean")

X=data.frame(sample_data(ps1))
x=data.frame(Body_site,Village,Disease,BMI,Subtribe,Shoes_type,Age,
Pets,Smoking,Gender,Antifungals,Presence_of_lesion,Family_history,Drugs,oilments)

G <- gower(D)
lambda <- eigen(G, only.values = T)$values
mdmr.res2 <- mdmr(X =x, G = G, lambda = lambda)
b=summary(mdmr.res2)

b=data.frame(b)
colnames(b)=c("statistic","DF","R2","p_val")
b$q_val=p.adjust(b$p_val,method="BH")
write.csv(b,"mdmr_all_1703_1000.csv")



b=read.csv("mdmr_oa_1703_1000.csv")
c=read.csv("adonis_oa_1703_1000.csv")

b=b[order(b$R2),]


p1 <- ggplot(data=b, aes(x=factor(Variable,levels=unique(b$Variable)), y=R2, fill=(q_val<0.05))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey")) +ggtitle("MDMR")+
  theme_classic() +labs(x="Variable",y="Effect sizes",fill="adjusted p-value")+
 coord_flip()+
  theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
plot.title=element_text(size=8,face="bold"),
                  legend.position = "bottom")




p2<-ggplot(data=c, aes(x=factor(Variable,levels=unique(b$Variable)), y=R2, fill=(q_val<0.05))) +
  geom_bar(stat="identity")+
scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))+
theme_classic()+ggtitle("Beta-diversity")+
 labs(x="Variable",y="Effect sizes",fill="adjusted p-value")+
  coord_flip()+
  theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
plot.title=element_text(size=8,face="bold"),
                  legend.position = "bottom")


q=ggarrange(p1+theme(legend.position="none",axis.title.x=element_blank()),p2+theme(legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank()),
nrow=1,ncol=2,widths=c(1.5,0.8),align="h")


c=read.csv("summary_interaction_oa_1703.csv")
c=filter(c, Variable != "(Omnibus)")
c=c[order(c$Pseudo.R2),]


p2<-ggplot(data=c, aes(x=factor(Variable,levels=c$Variable), y=Pseudo.R2, fill=(Analytic.p.value<0.05))) +
  geom_bar(stat="identity")+
scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))+
theme_classic()+ggtitle("MDMR  mixed models")+
 labs(x="Variable",y="Effect sizes",fill="p-value")+
  coord_flip()+
  theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
plot.title=element_text(size=8,face="bold"),
                  legend.position = "right")+ylim(0,0.05)



q2=ggarrange(q,p2+theme(legend.position="none"),ncol=1,nrow=2,heights=c(1,0.4))



#all


b=read.csv("mdmr_all_1703_1000.csv")
c=read.csv("adonis_all_1703_1000.csv")

b=b[order(b$R2),]


p1 <- ggplot(data=b, aes(x=factor(Variable,levels=unique(b$Variable)), y=R2, fill=(q_val<0.05))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey")) +ggtitle("MDMR")+
  theme_classic() +labs(x="Variable",y="Effect sizes",fill="adjusted p-value")+
 coord_flip()+
  theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
plot.title=element_text(size=8,face="bold"),
                  legend.position = "bottom")




p2<-ggplot(data=c, aes(x=factor(Variable,levels=unique(b$Variable)), y=R2, fill=(q_val<0.05))) +
  geom_bar(stat="identity")+
scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))+
theme_classic()+ggtitle("Beta-diversity")+
 labs(x="Variable",y="Effect sizes",fill="adjusted p-value")+
  coord_flip()+
  theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
plot.title=element_text(size=8,face="bold"),
                  legend.position = "bottom")


q3=ggarrange(p1+theme(legend.position="none",axis.title.x=element_blank()),p2+theme(legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank()),
nrow=1,ncol=2,widths=c(1.5,0.8),align="h")



#all
c=read.csv("summary_interaction_1703_all.csv")
c=filter(c, Variable != "(Omnibus)")
c=c[order(c$Pseudo.R2),]


p2<-ggplot(data=c, aes(x=factor(Variable,levels=c$Variable), y=Pseudo.R2, fill=(Analytic.p.value<0.05))) +
  geom_bar(stat="identity")+
scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"))+
theme_classic()+ggtitle("MDMR  mixed models")+
 labs(x="Variable",y="Effect sizes",fill="p-value")+
  coord_flip()+
  theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
plot.title=element_text(size=8,face="bold"),
                  legend.position = "right")+ylim(0,0.05)


q4=ggarrange(q3,p2+theme(legend.position="none"),ncol=1,nrow=2,heights=c(1,0.4))



q5=ggarrange(q2,q4,nrow=1,ncol=2)


ggsave("Fig_S3A-B.svg",q5,height=7,width=7.7)


