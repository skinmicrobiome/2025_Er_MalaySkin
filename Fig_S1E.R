


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



df=df %>% select(4,6:8,10:11)

n2=melt(df)


n2$value=as.factor(n2$value)


n2$variable=factor(n2$variable,levels=c("Miconazole","Ketoconazole",
"Fluconazole","Itraconazole","Posaconazole","Voriconazole","Isavuconazole",
"Terbinafine","Amphotericin.B","Griseofulvin","Ciclopirox","Micafungin",
"Olorofim"))


n2$variable=recode(n2$variable,"Miconazole"="MCZ","Ketoconazole"="KTZ",
"Fluconazole"="FLZ","Itraconazole"="ITC","Posaconazole"="PCZ",
"Voriconazole"="VOR","Isavuconazole"="ISA","Amphotericin.B"="AMB",
"Griseofulvin"="GRI","Ciclopirox"="CPX","Micafungin"="MIC","Olorofim"="OLO")




n2$value=recode(n2$value,"2"="2.000","1"="1.000","0.5"="0.500","4"="4.000","16"="16.000",
"0.25"="0.250","8"="8.000")



#scale_y_continuous(labels=scales::percent) +
myplot <- ggplot(n2, aes(value,fill=value)) + 
          geom_bar(aes(y = (..count..)/sum(..count..)),fill="#9999ff") +
       geom_line(stat="count", aes(y = (..count..)/sum(..count..)), color="red", group=1,linewidth=1) +          
xlab("MIC value/ μg/ml")+
  ylab("Percentage")+ theme_classic()+scale_fill_jco()+theme(legend.position="none")+
facet_grid(~variable,scale="free",space="free")


myplot <- ggplot(n2, aes(value, fill=value)) + 
          geom_bar(fill="#9999ff") +
          geom_line(stat="count", color="red", group=1, linewidth=1) +          
          xlab("MIC value (μg/ml)") +
          ylab("Isolates") + 
          theme_classic() +
                  theme(legend.position="none") +
          facet_grid(~variable, scale="free", space="free")

a=myplot+
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
        axis.title.x = element_text(size=8),
strip.text.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none",plot.title = element_text(size=22))+ guides(colour = guide_legend(override.aes = list(size=5)))


ggsave("Fig_S1E.svg",a,height=1.5,width=8.1)


