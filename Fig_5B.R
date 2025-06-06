
library(readxl)
library(knitr)

library(readxl)
library(knitr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(XML)
library(stringr)
library(phyloseq)
library(vegan)
library(microbiome)


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

bracken=readRDS("relab_cory_recoded.rda")

ps=subset_samples(bracken, site_specific %in% c("Tw"))
df=psmelt(ps)

df2=df[order(-df$Abundance),]

x=head(unique(df2$ta5),n=19)


df2 <- df2 %>%
  mutate(sp3 = ifelse(ta5 %in% x , ta5, "others"))


library(dplyr)


df2 <- df2 %>%
  mutate(sp3 = factor(sp3, levels = c("others", setdiff(unique(ta5[ta5 %in% x]), "ta5"))))


df2$village=factor(df2$village,levels=c("Kuala Lumpur","Washington DC","Village F",
                                        "Village D","Village G","Village AB",
                                        "Village C"))


library(dplyr)


df3=filter(df2, Abundance >0)

df3$sp3=factor(df3$sp3,levels=c("Novel Malaysia",
"Corynebacterium tuberculostearicum",
"Corynebacterium kefirresidentii",
"Corynebacterium hallux"  ,   "SMGC",                                   
"Corynebacterium matruchotii ATCC 14266",                       
"Corynebacterium ureicelerivorans",                     
"Corynebacterium nuruki S6-4",        
"Corynebacterium striatum",                   
"Corynebacterium propinquum", 
"Corynebacterium massiliense DSM 45435",
"Corynebacterium macclintockiae",
"Corynebacterium durum",
"Corynebacterium sp900232865",
"Corynebacterium lujinxingii",
"Corynebacterium amycolatum",
"Corynebacterium faecigallinarum",
"Corynebacterium incognita",
"Corynebacterium otitidis ATCC 51513",                                                
"Corynebacterium pseudodiphtheriticum DSM 44287", 
"Corynebacterium jeddahense",
"Corynebacterium ihumii",
"Corynebacterium mycetoides",
"Corynebacterium accolens",
"Corynebacterium resistens DSM 45100",
"Corynebacterium appendicis CIP 107643",
"Corynebacterium sanguinis",
"Corynebacterium efficiens YS-314",
"Corynebacterium variabile",
"others" ,                                 
"Novel OA"                         
)  )



colors <- c("darkolivegreen","peru", "tan4", "sandybrown","sienna3", "#E6E6E6", "#DCDCDC", "#F5F5F5", "#FAFAFA",
"#F8F8FF", "#F0F8FF", "#F0FFFF", "#E0FFFF", "#E6F3FF",
"#F0F3FF", "#FFF0F5", "#FBEDF2", "#BEC0CB", "#FFF5F5",
"#FFF0F0", "#FFE4E1", "#FFEBF0", "#FFE6EE", "#FFF0FF",
"#F5E0FF", "#F0E6FF", "#E6E6FA", "#F8F0FF", "#EAEAEA", "palegreen","seagreen")


taxa_colors <- setNames(colors, levels(df3$sp3))

# Print the result
print(taxa_colors)


d=ggplot(df3[df3$site_specific == 'Tw', ]) +
  aes(x=subject_id, y=Abundance, fill=sp3) +
  geom_bar(stat="identity", position="fill", width=0.8) +
  theme_bw(base_size=12) +ggtitle("Toe web")+
  ylab('Relative Abundance/%') +xlab("Samples")+
  #scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0', '', '50', '', '100'), limits=c(0,1)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=8,face="bold"),
        axis.title.y=element_text(face="bold",size=8), axis.text.y=element_text(face="bold",size=8)) +
  theme(legend.key.size=unit(0.5, "cm"), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(legend.position='right') +
  theme(legend.box="vertical", legend.title=element_text(), legend.text=element_text()) +
  scale_fill_manual("Classification", values=taxa_colors) +
  theme(strip.text.x=element_text(), strip.background=element_rect(fill="white")) +
facet_grid(~village,scale="free",space='free')
  

d2=d+theme(axis.text.x = element_blank(),
           axis.text.y = element_text(size=8,face="bold"),
           axis.title.x = element_text(size=8,face="bold"),
           axis.title.y = element_text(size=8,face="bold"),
           strip.text.x = element_blank(),
           legend.text = element_text(size=8,face="italic"),
           legend.title = element_text(size=8))

c=d2+theme(strip.text.x =element_text(size=8,face="bold"),strip.text.y = element_text(
    size = 8, face="bold",angle=270))


c=c+theme(legend.position="right")+  guides(fill = guide_legend(reverse=FALSE, nrow=30))

ggsave("Fig_5B.svg",c,height=9,width=13)




