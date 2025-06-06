


# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)

library(readxl)
library(knitr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(XML)
library(RcppParallel)
library(stringr)
library(phyloseq)
library(vegan)
library(rstatix)
library(ggpubr)
library(pacman)
library(reshape2)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


df=read.csv("resistance_150125.csv")

df=df %>% select(2,4)


n2=melt(df)


n2$value=as.factor(n2$value)


n2$value=recode(n2$value,"2"="2.000","1"="1.000","0.5"="0.500","4"="4.000")


n2$variable=factor(n2$variable, levels=c("Terbinafine","Griseofulvin"))


n2$variable=recode(n2$variable,"Fluconazole"="Flz")

myplot <- ggplot(n2, aes(value, fill=value)) + 
          geom_bar(fill="#9999ff") +
          geom_line(stat="count", color="red", group=1, linewidth=1) +          
          xlab("MIC value/ μg/ml") +
          ylab("Number of isolates") + 
          theme_classic() +
                  theme(legend.position="none") +
          facet_wrap(~variable, scale="free_x",ncol=3)


a=myplot+
  theme(axis.text.y = element_text(size=8,face="bold"),
        axis.text.x = element_text(size=8,face="bold",angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_text(size=8,face="bold"),
strip.text.x = element_text(size=8,face="bold"),
        axis.title.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))


ggsave("Fig_2C_2004.jpg",a,height=9,width=10)


