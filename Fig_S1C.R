
setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

library(pacman)
library(viridis)

pacman::p_load(dplyr,phyloseq,tidyverse, microbiome,treeio,ggtree,ape,ggnewscale)


df2=read.csv("parsnp_2702.csv")

df2$IND=ifelse(df2$IND_1 == df2$IND_2,"Within individual","no")

#df2=filter(df, Distance >0)

df2$VILL_1 <- gsub("\\d+", "", df2$IND_1)
df2$VILL_2 <- gsub("\\d+", "", df2$IND_2)


df2$within_village <- ifelse(df2$Family == "ND" & 
                           df2$IND == "no" & 
                           df2$VILL_1 == df2$VILL_2, 
                           "Within village", 
                           "no")

df2$between_village=ifelse(df2$VILL_1 == df2$VILL_2,"no","Between village")



df2$within_family <- ifelse(df2$Family %in% c("C_002", "C_003", "F_001") & 
                         df2$IND == "no", 
                         "within family",
                         "no")

df2$summary <- ifelse(df2$IND == "Within individual", "Within individual",
                    ifelse(df2$within_village == "Within village", "Within village",
                          ifelse(df2$between_village == "Between village", "Between village",
                                ifelse(df2$within_family == "within family", "Within family",
                                      "no"))))



df2$Family_1 <- case_when(
  df2$IND_1 %in% c("F09", "F16") ~ "Family 2",
  df2$IND_1 %in% c("C07", "C08") ~ "Family 3",
  df2$IND_1 %in% c("C05", "C03", "C39", "C40") ~ "Family 1",
  TRUE ~ "ND"
)



df2$Family_2 <- case_when(
  df2$IND_2 %in% c("F09", "F16") ~ "Family 2",
  df2$IND_2 %in% c("C07", "C08") ~ "Family 3",
  df2$IND_2 %in% c("C05", "C03", "C39", "C40") ~ "Family 1",
  TRUE ~ "ND"
)




df2$Family_1=factor(df2$Family_1, levels=c("Family 1","Family 2","Family 3","ND"))


df2$Family_2 = factor(df2$Family_2, levels=c("Family 1", "Family 2", "Family 3", "ND"))

df2$blank="blank"
# For p2, you need to make sure it uses the same Genome1 ordering as in p1
p2 <- ggplot(df2, aes(Genome1, blank, fill=Family_1)) +
  geom_tile() + theme_classic() +
  xlab("") + ylab("") + labs(fill="Family") +
  scale_fill_manual(values=c("red", "blue", "purple", "honeydew")) +
  theme(axis.text.x = element_text(size=7, angle=90, vjust=0.5, hjust=0.5),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank())


p1<-ggplot(df2,aes(Genome1,Genome2,fill=Distance)) +
  geom_tile() +
  xlab("") + ylab("") +
  scale_fill_viridis(option="turbo",discrete=FALSE,limits = c(0, 500)) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
axis.title.x = element_blank(),axis.title.y = element_blank(),axis.ticks = element_blank())




p3<-ggplot(df2,aes(blank,Genome2,fill=Family_2)) +
  geom_tile() +theme_classic()+
  xlab("") + ylab("") +
  scale_fill_manual(values=c("red","blue","purple","honeydew")) +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size=6),
axis.title.x = element_blank(),axis.title.y = element_blank())

q=ggarrange(p3+theme(legend.position="none"),p1+theme(legend.position="none"), widths=c(0.3,1),nrow=1,ncol=2)

q2=ggarrange(NULL,p2+theme(legend.position="none"),widths=c(0.3,1),nrow=1,ncol=2)


q3=ggarrange(q,q2,nrow=2,ncol=1,heights=c(1,0.3),align="h")


p1=p1+theme(legend.position="right")

p2=p2+theme(legend.position="right")

p3=p2+theme(legend.position="bottom")

q4=ggarrange(get_legend(p1),get_legend(p2),widths=c(0.4,1),nrow=1,ncol=2,
align="hv")

ggsave("Fig_S1C_2702.svg",q3,height=4,width=5)







