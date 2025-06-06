
setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

library(pacman)


pacman::p_load(dplyr,phyloseq,tidyverse, microbiome,treeio,ggpubr,reshape2,ggnewscale)

data2=read.csv("resist_0803.csv")

data2=filter(data2, remark == "ok")



df=data2 %>% select(6,7:12,13,14,15,16)
write.csv(df, "resist_coded.csv")

z1=melt(df)


z1$variable=factor(z1$variable,levels=c("Miconazole","Ketoconazole","Fluconazole",
"Itraconazole","Voriconazole","Griseofulvin","Terbinafine","Ciclopirox"))


z2=filter(z1, variable =="Miconazole")

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,4),labels=c("0","4"),
                       limits=c(0,4))

a<- x+
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
              legend.position = "none")




z2=filter(z1, variable =="Ketoconazole")

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,8),labels=c("0","8"),
                       limits=c(0,8))

b<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")


z2=filter(z1, variable =="Fluconazole")

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,16),labels=c("0","16"),
                       limits=c(0,16))

c<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")



z2=filter(z1, variable =="Itraconazole")


summary(z2$value)

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,2),labels=c("0","2"),
                       limits=c(0,2))

d<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")



z2=filter(z1, variable =="Voriconazole")


summary(z2$value)

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,1),labels=c("0","1"),
                       limits=c(0,1))

e<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")



z2=filter(z1, variable =="Griseofulvin")


summary(z2$value)

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,4),labels=c("0","4"),
                       limits=c(0,4))

f<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")




z2=filter(z1, variable =="Terbinafine")


summary(z2$value)

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,2),labels=c("0","2"),
                       limits=c(0,2))

g<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")



z2=filter(z1, variable =="Ciclopirox")


summary(z2$value)

x<-ggplot(z2, aes(x=variable,y = NEW_ID)) + 
  geom_tile(aes(fill= value)) +theme_classic()+
  labs(x="Isolates", y="Antifungal", 
       fill="MIC values")+geom_text(aes(label=sprintf("%.2f", value)),size=2.3) +
  theme(axis.text.y=element_text(size=8), 
        legend.position = "bottom")+
  scale_fill_gradientn(colours = c("white", "red"),
                       breaks=c(0,1),labels=c("0","1"),
                       limits=c(0,1))

h<- x+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
              legend.position = "none")

q=ggarrange(a,b,c,d,e,f,g,h,nrow=1,ncol=8,widths=c(1.2,0.7,0.7,0.7,0.7,0.7,0.7,0.7),align="h")

a1 <- a + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
b1 <- b + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
c1 <- c + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
d1 <- d + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
e1 <- e + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
f1 <- f + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
g1 <- g + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
h1 <- h + theme(plot.margin = margin(0, 0, 0, 0, "pt"))



q = ggarrange(a1, b1, c1, d1, e1, f1, g1, h1,
              nrow = 1,
              ncol = 8,
              widths = c(1.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
              align = "h"
               )  # Reduce padding to minimum



ggsave("resist_0803.svg",q,height=7.6,width=4)


