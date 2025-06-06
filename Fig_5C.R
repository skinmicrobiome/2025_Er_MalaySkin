
library(readxl)
library(knitr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(microbiome)
library(stringr)
library(phyloseq)
library(vegan)
library(pacman)

pacman::p_load( MicrobiomeStat,tidyverse, dplyr, ggplot2,ggsci,vegan,reshape2,phyloseq,ggpubr,readr,tidyverse)

setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


a=read.csv("gene_data_cory_1506_coded.csv")

info2=read.csv("Cory_info_2805_coded.csv")
a2=filter(a, !description %in% c("hypothetical protein","putative protein"))
a3 <- a2[!grepl("location:", a2$description, ignore.case = TRUE), ]

info2=read.csv("Cory_info_2805_coded.csv")

a3$id2=a3$gff_file
a4=left_join(a3,info2)

a5=filter(a4, id2 %in% info2$id2)


n=data.frame(table(a5$id2,a5$description))


# Correcting the approach to create a dataframe with n$Var1 as colnames and n$Var2 as rownames
n_df <- dcast(n, Var2 ~ Var1, value.var = "Freq")
rownames(n_df) <- n_df$Var2
n_df <- n_df[,-1] # Remove the Var2 column after setting it as rownames

data4=filter(info2, location != "Reference")
data5=filter(data4, id2 %in% colnames(n_df))

data5$oauc=recode(data5$location,"Kuala_Lumpur"="Urban","Washington_DC"="Urban",.default="Orang Asli")

rownames(data5)=data5$id2

linda.obj  <- linda(n_df, data5, formula = '~oauc', corr.cut=0,
           feature.dat.type = 'count',  pseudo.cnt = 0.5,          prev.filter = 0.05, is.winsor = TRUE, outlier.pct = 0.03,
           p.adj.method = "BH", alpha = 0.05)


linda.plot(linda.obj, c('oaucUrban'),
           titles = c('OAUC: Urban v.s. Orang Asli'), 
           alpha = 0.1, lfc.cut = 1, legend = TRUE, directory = NULL,
            width = 11, height = 8)

c=linda.obj$output 
c=data.frame(c)
colnames(c)=c("baseMean","log2FoldChange","lfcSE","stat",
"p_val","q_val","reject","df")
c$pathway=rownames(c)

# Create a new column named "Category"
c$Category <- ifelse(c$log2FoldChange > 0.5 & c$p_val < 0.05, "Urban",
                      ifelse(c$log2FoldChange < -0.5 & c$p_val < 0.05, "Orang Asli", "No significant"))


c$Category=factor(c$Category, levels=c("Orang Asli","Urban","No significant"))

write.csv(c,"linda_cory_0507_2.csv")

c=read.csv("linda_cory_0507_2.csv")


c$Category=factor(c$Category, levels=c("Orang Asli","Urban","No significant"))

c$adjust_p=p.adjust(c$p_val, method="BH")


c2=filter(c, q_val <0.05)


table(c$Category)

p1 <- c %>%
  ggplot(aes(x = log2FoldChange, y = -log10(p_val), color = Category)) + 
  geom_point() +
  scale_color_manual(values = c("seagreen", "sienna2", "grey")) +
  theme_classic() +
  labs(x = "bias-corrected coefficients", y = "-log10 (p-value)", color = "Cohort") +
  geom_text(data = subset(c, Category %in% c("Orang Asli", "Urban")), aes(label = pathway), vjust = 1.5, hjust = 0.5, size = 2, check_overlap = TRUE)+
 theme( axis.text.x = element_text(size=8,angle = 90, vjust = 0.5,face="bold"),
         axis.text.y = element_text(size=8,face="bold"),
         axis.title.x = element_text(size=8,face="bold"),
         axis.title.y = element_text(size=8,face="bold"),
         legend.text = element_text(size=8),
         legend.title = element_text(size=8),
         legend.position = "bottom",plot.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=5)))
p1


ggsave("Fig_5C.svg",p1+theme(legend.position="none"),height=3,width=4)




