
library(dplyr)
library(nloptr)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(readr)


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

ps=readRDS("phyloseq_0505.rda")

sample_data(ps)$oauc=recode(sample_data(ps)$village,"Washington DC"="UC",
                               "Kuala Lumpur"="UC",
                             .default="OA")

ps1=subset_samples(ps, oauc == "UC")

otu4=data.frame(otu_table(ps1))


otu4$Time <- apply(otu4, 1, function(row) {
  ifelse(any(row > 0), which(row > 0)[1], NA)
})

n=data.frame(otu4$Time)
tax=data.frame(tax_table(ps1))
tax$sp2=paste(tax$Genus,tax$Species)
n$taxa=tax$sp2
n$number=rownames(n)
colnames(n)=c("Time","Taxa","Number")

library(ggplot2)

n$Time=as.factor(n$Time)

n3=data.frame(table(n$Time))
n3$Cum <- cumsum(n3$Freq)

n3$Cum=as.numeric(as.character(n3$Cum))
n3$Var1=as.numeric(as.character(n3$Var1))

n3$Group="UC"

ps2=subset_samples(ps, oauc == "OA")

otu4=data.frame(otu_table(ps2))


otu4$Time <- apply(otu4, 1, function(row) {
  ifelse(any(row > 0), which(row > 0)[1], NA)
})

n=data.frame(otu4$Time)
tax=data.frame(tax_table(ps2))
tax$sp2=paste(tax$Genus,tax$Species)
n$taxa=tax$sp2
n$number=rownames(n)
colnames(n)=c("Time","Taxa","Number")

library(ggplot2)

n$Time=as.factor(n$Time)

n4=data.frame(table(n$Time))
n4$Cum <- cumsum(n4$Freq)

n4$Cum=as.numeric(as.character(n4$Cum))
n4$Var1=as.numeric(as.character(n4$Var1))

n4$Group="OA"



df=rbind(n3,n4)


p3 <- ggplot(df, aes(x = Var1, y = Cum,color=Group)) +
 # geom_line(linewidth=0) +
  theme_classic() +
  labs(x = "Number of samples", y = "Number of Species",color="Group")+
geom_smooth(linewidth=1,method="loess")+
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 1000))+
  scale_color_manual(values=c("seagreen","sienna2"))+
  theme( axis.text.x = element_text(size=8),
         axis.text.y = element_text(size=8),
         axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         legend.text = element_text(size=8),
         legend.title = element_text(size=8),
         legend.position = "bottom")+
  guides(colour = guide_legend(nrow=1,override.aes = list(size=5)))

ggsave("Fig_S3C.svg",p3,dpi=320,height=3,width=4)



