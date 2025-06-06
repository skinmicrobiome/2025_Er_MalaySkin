
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


setwd("/Users/ery2/Desktop/QA_TIDY_PRE20/smgc_malaysia/ref_smgc/phylo_form")


metadata_finder = function(df, metadata){
  subject_id = NULL
  subject_de = NULL
  group_id = NULL
  site = NULL
  for (i in as.character(df$variable)){
    subject_id = append(subject_id, metadata$subject_id[metadata$library_id == i])
    group_id <- append(group_id, metadata$group_id[metadata$library_id == i])
    site = append(site, metadata$site_specific[metadata$library_id== i])
  }
  df$subject_id <- subject_id
  df$cohort <- group_id
  df$site <- site
  return(df)
}




ps=readRDS("phyloseq_0505.rda")

bracken.phy=subset_samples(ps,!village %in% c("Washington DC","Kuala Lumpur"))

bracken.phy=subset_samples(bracken.phy, site_specific != "Ctrl")
bracken.phy <- subset_taxa(bracken.phy, !Species %in% c('sapiens', 'gondii', 'mexicana', 'bigemina'))

prev <- apply(X=otu_table(bracken.phy), MARGIN=ifelse(taxa_are_rows(bracken.phy), yes=1, no=2), FUN=function(x){sum(x > 0)})
df_prev <- data.frame(Prevalence=prev, TotalAbundance=taxa_sums(bracken.phy), tax_table(bracken.phy))
low_phylum <- plyr::ddply(df_prev, 'Phylum', function(df1){cbind(mean(df1$Prevalence), sum(df1$Prevalence), mean(df1$TotalAbundance), sum(df1$TotalAbundance))})
filterPhyla <- low_phylum$Phylum[low_phylum[2] <= 0.01*nsamples(bracken.phy)]
bracken.phy <- subset_taxa(bracken.phy, !Phylum %in% filterPhyla)

# convert OTU table to relative abundances
ra <- transform_sample_counts(bracken.phy, function(x) x/sum(x))
df <- cbind(data.frame(tax_table(ra)), data.frame(otu_table(ra)))

# sort the table by total relative abundances of species
df_sorted <- df[order(rowSums(df[, 8:ncol(df)]), decreasing=TRUE), ]

ccca_sorted=df_sorted
ccca_meta=sample_data(bracken.phy)



TargetLvl <- function(df, mytab){
  joinby <- 'library_id'
  target_taxa <- c('Bacteria (Others)',
                   'Actinomycetota',
                    'Bacillota',
                             'Pseudomonadota',
                  'Eukaryota (Others)','Malassezia','Trichophyton',
                   'Viruses')
  df_origin <- df
  df$Species <- gsub('_', ' ', df$Species)
  taxa_only <- df[,c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
  df <- melt(df, id=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
  df$Target <- df$Kingdom
  df$Target[df$Kingdom == 'Archaea'] <- 'Bacteria (Others)' 
  df$Target[df$Kingdom == 'Bacteria'] <- 'Bacteria (Others)'
  df$Target[df$Phylum == 'Actinomycetota'] <- 'Actinomycetota'
  
  
  df$Target[df$Phylum == 'Bacillota'] <- 'Bacillota'

  df$Target[df$Phylum == 'Pseudomonadota'] <- 'Pseudomonadota'

  df$Target[df$Kingdom == 'Eukaryota'] <- 'Eukaryota (Others)'
  df$Target[df$Genus == 'Malassezia'] <- 'Malassezia'
  df$Target[df$Genus == 'Trichophyton'] <- 'Trichophyton'

  
  df$Target[df$Kingdom == 'Viruses'] <- 'Viruses'
  
  
  
  
  
  df$Target <- factor(df$Target, levels=target_taxa)
  df$Kingdom <- NULL
  df$Phylum <- NULL
  df$Class <- NULL
  df$Order <- NULL
  df$Family <- NULL
  df$Genus <- NULL
  df$Species <- NULL
  
  df <- data.frame(df %>% group_by(variable, Target) %>% summarise_all(list(sum)))
  df <- metadata_finder(df, mytab)
  return(df)
}


df <- TargetLvl(ccca_sorted, ccca_meta)


color_pall <- c(brewer.pal(9, 'YlOrRd')[c(1)],
                brewer.pal(9, 'Greens')[c(7)],
                brewer.pal(9, 'Blues')[c(3)],
                brewer.pal(9, 'Oranges')[c(6)],
                brewer.pal(9, 'Greys')[c(3,5,7)],
                brewer.pal(9, 'PuRd')[c(8)])


df$site=factor(df$site,levels=c("Fh","Vf","Tw","Ctrl"))
levels(df$site)=c("Foreheads","Volar forearm","Toewebs","Control")


n=data.frame(ccca_meta$subject_id,ccca_meta$library_id,ccca_meta$village,ccca_meta$Disease,
             ccca_meta$Group, ccca_meta$a_u)
colnames(n)=c("id","variable","village","disease","group","status")

df2=dplyr::left_join(df,n)

df2$disease=recode(df2$disease, "Negative"="Healthy control")


df4=filter(df2, disease %in% c("Scabies"))


mean_values <- aggregate(df4$value, by=list(df4$Target), FUN=mean)
colnames(mean_values) <- c("Target", "MeanValue")

data1=data.frame(mean_values)
data1$Disease="Scabies"


df4=filter(df2, disease %in% c("Healthy control"))


mean_values <- aggregate(df4$value, by=list(df4$Target), FUN=mean)
colnames(mean_values) <- c("Target", "MeanValue")

data2=data.frame(mean_values)
data2$Disease="Healthy control"


df_3=rbind(data1,data2)


df_3$Target=factor(df_3$Target,levels=c("Bacteria (Others)","Actinomycetota","Bacillota","Pseudomonadota",
                                        "Viruses","Eukaryota (Others)","Malassezia","Trichophyton"))

a=ggplot(df_3, aes(x = Disease, y = MeanValue*100, fill = Target)) +
  geom_col() +
  scale_fill_manual(values=c("yellow3", "bisque", "steelblue", "seagreen", "grey2", "salmon2", "thistle2", "violetred2")) +
  coord_polar("y")+theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=0.5,vjust=0.5,size=8,face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom")+ guides(colour = guide_legend(override.aes = list(size=5)))




#h tv

df4=filter(df2, disease %in% c("Tinea versicolor"))


mean_values <- aggregate(df4$value, by=list(df4$Target), FUN=mean)
colnames(mean_values) <- c("Target", "MeanValue")

data1=data.frame(mean_values)
data1$Disease="Tinea versicolor"


df4=filter(df2, disease %in% c("Healthy control"))


mean_values <- aggregate(df4$value, by=list(df4$Target), FUN=mean)
colnames(mean_values) <- c("Target", "MeanValue")

data2=data.frame(mean_values)
data2$Disease="Healthy control"


df_3=rbind(data1,data2)


df_3$Target=factor(df_3$Target,levels=c("Bacteria (Others)","Actinomycetota","Bacillota","Pseudomonadota",
                                        "Viruses","Eukaryota (Others)","Malassezia","Trichophyton"))

b=ggplot(df_3, aes(x = Disease, y = MeanValue*100, fill = Target)) +
  geom_col() +
  scale_fill_manual(values=c("yellow3", "bisque", "steelblue", "seagreen", "grey2", "salmon2", "thistle2", "violetred2")) +
  coord_polar("y")+theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=0.5,vjust=0.5,size=8,face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom")+ guides(colour = guide_legend(override.aes = list(size=5)))


#ti

df4=filter(df2, disease %in% c("Tinea imbricata"))


mean_values <- aggregate(df4$value, by=list(df4$Target), FUN=mean)
colnames(mean_values) <- c("Target", "MeanValue")

data1=data.frame(mean_values)
data1$Disease="Tinea imbricata"


df4=filter(df2, disease %in% c("Healthy control"))


mean_values <- aggregate(df4$value, by=list(df4$Target), FUN=mean)
colnames(mean_values) <- c("Target", "MeanValue")

data2=data.frame(mean_values)
data2$Disease="Healthy control"


df_3=rbind(data1,data2)


df_3$Target=factor(df_3$Target,levels=c("Bacteria (Others)","Actinomycetota","Bacillota","Pseudomonadota",
                                        "Viruses","Eukaryota (Others)","Malassezia","Trichophyton"))

c=ggplot(df_3, aes(x = Disease, y = MeanValue*100, fill = Target)) +
  geom_col() +
  scale_fill_manual(values=c("yellow3", "bisque", "steelblue", "seagreen", "grey2", "salmon2", "thistle2", "violetred2")) +
  coord_polar("y")+theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=0.5,vjust=0.5,size=8,face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8,face="bold"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom")+   
guides(fill = guide_legend(ncol=2,override.aes = list(size=5)))




q2=ggarrange(c+theme(legend.position="none"),b+theme(legend.position="none"),
            a+theme(legend.position="none"),nrow=3,ncol=1,heights=c(1,1,1))



ggsave("Fig_6B.svg",q2, height=9, width=4)

ggsave("Fig_6B_legend.svg",get_legend(c), height=3, width=4)



