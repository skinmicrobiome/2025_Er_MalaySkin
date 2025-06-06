
setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

library(pacman)


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




df2$summary=factor(df2$summary,levels=c("Between village","Within village","Within family","Within individual"
))


df21=filter(df2, Distance >0)
p5<-ggplot(df21,aes(Distance)) +
    xlab("Pairwise SNP distance (ParSNP)") +theme_classic()+
ylab("Number of comparisons")+
  geom_histogram(alpha=1,fill="#cf3379") +
  facet_grid(summary~ blank,scales="free_y")+
 theme( axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
         axis.text.y = element_text(size=8),
axis.title.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         legend.position="none",strip.text.x=element_blank()
       )+ylim(0,200)


p5

ggsave("Fig_S1B.svg",p5,height=5,width=4)









