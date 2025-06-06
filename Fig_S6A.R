
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


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

#bracken.phy <- import_biom('kraken_arc.biom', parseFunction=parse_taxonomy_greengenes)


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


# Create sorted relative abundance table from .biom file

# import your biom file



ps=readRDS("phyloseq_0505.rda")



ps=subset_samples(ps,village != "Washington DC")

bracken.phy=subset_samples(ps,village != "Kuala Lumpur")
#bracken.phy=subset_samples(bracken.phy, Disease != "Negative")
bracken.phy=subset_samples(bracken.phy, site_specific != "Ctrl")

bracken.phy <- subset_taxa(bracken.phy, !is.na(Phylum))
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
  target_taxa <- c(
                   'Actinomycetota (Others)',
                   'Cutibacterium','Corynebacterium', 'Pseudomonadota (Others)',
                   'Moraxella', 'Acinetobacter', 
                   'Paracoccus','Pseudomonas','Bacillota (Others)',
                   'Staphylococcus','Streptococcus','Eukaryota (Others)','Trichophyton concentricum','Malassezia','M. japonica',
                   'Others')
  df_origin <- df
  df$Species <- gsub('_', ' ', df$Species)
  taxa_only <- df[,c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
  df <- melt(df, id=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
  df$Target <- df$Kingdom
  df$Target[df$Kingdom == 'Archaea'] <- 'Others' 
  df$Target[df$Kingdom == 'Bacteria'] <- 'Others'
  df$Target[df$Phylum == 'Actinomycetota'] <- 'Actinomycetota (Others)'
  df$Target[df$Genus == 'Cutibacterium'] <- 'Cutibacterium'
  df$Target[df$Genus == 'Corynebacterium'] <- 'Corynebacterium'
#  df$Target[df$Genus == 'Corynebacterium'& df$Species == 'resistens' ] <- 'C. resistens'

  
  df$Target[df$Phylum == 'Bacillota'] <- 'Bacillota (Others)'
  df$Target[df$Genus == 'Staphylococcus'] <- 'Staphylococcus'
  df$Target[df$Genus == 'Streptococcus'] <- 'Streptococcus'

  df$Target[df$Phylum == 'Pseudomonadota'] <- 'Pseudomonadota (Others)'
  df$Target[df$Genus == 'Moraxella'] <- 'Moraxella'
  df$Target[df$Genus == 'Acinetobacter'] <- 'Acinetobacter'
  df$Target[df$Genus == 'Paracoccus'] <- 'Paracoccus'
  df$Target[df$Genus == 'Pseudomonas'] <- 'Pseudomonas'
  
  df$Target[df$Kingdom == 'Eukaryota'] <- 'Eukaryota (Others)'
  df$Target[df$Genus == 'Malassezia'] <- 'Malassezia'
  df$Target[df$Genus == 'Malassezia' & df$Species == 'japonica'] <- 'M. japonica'
  df$Target[df$Genus == 'Trichophyton'] <- 'Trichophyton concentricum'

  
  df$Target[df$Kingdom == 'Viruses'] <- 'Others'
  
  
  
  
  
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



color_pall <- c(
                 brewer.pal(9, 'Oranges')[c(1,4,6)],      
                brewer.pal(9, 'BuGn')[c(9,8,7,6,5)],
                brewer.pal(9, 'Blues')[c(3,9,7)],
                brewer.pal(9, 'PuRd')[c(1,7,3,5)],
 brewer.pal(9, 'YlOrRd')[c(3)]
               )



df$site=factor(df$site,levels=c("Fh","Vf","Tw","Ctrl"))
levels(df$site)=c("Forehead","Volar forearm","Toeweb","Control")


n=data.frame(ccca_meta$subject_id,ccca_meta$library_id,ccca_meta$village,ccca_meta$Disease,
             ccca_meta$Group, ccca_meta$a_u)
colnames(n)=c("id","variable","village","disease","group","status")

df2=dplyr::left_join(df,n)



df2$disease=recode(df2$disease, "Negative"="Healthy control")
df2$status=factor(df2$status,levels=c("A","U"))
levels(df2$status)=c("Nonlesional skin","Lesional skin")

df2=filter(df2, disease %in% c("Healthy control","Tinea imbricata"))

df2$group_status <- ifelse(df2$disease == "Healthy control", "Healthy control", ifelse(df2$disease == "Tinea imbricata", df2$status, NA))


df2$group_status=recode(df2$group_status,"1"="Lesional skin","2"="Nonlesional skin")

df2$group_status=factor(df2$group_status,levels=c("Healthy control","Nonlesional skin","Lesional skin"))


d=ggplot(df2[df2$site == 'Forehead', ]) +
  aes(x=subject_id, y=value, fill=Target,color=Target) +
  geom_bar(stat="identity", position="fill", width=0.8) +
  theme_bw(base_size=12) +
  ylab('Relative Abundance/%') +xlab("Status")+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0', '', '50', '', '100'), limits=c(0,1)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=8,face="bold"),
        axis.title.y=element_text(face="bold",size=8), axis.text.y=element_text(face="bold",size=8)) +
  theme(legend.key.size=unit(0.5, "cm"), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(legend.position='right') +
  theme(legend.box="vertical", legend.title=element_text(), legend.text=element_text()) +
  scale_fill_manual("Classification", values=color_pall) +scale_color_manual(values=color_pall)+
  theme(strip.text.x=element_text(), strip.text.y=element_text(angle=0, ), strip.background=element_rect(fill="white")) +
  guides(fill = guide_legend(reverse=FALSE, ncol=1))


d2=d+theme(axis.text.x = element_blank(),
           axis.text.y = element_text(size=8,face="bold"),
           axis.title.x = element_text(size=8,face="bold"),
           axis.title.y = element_text(size=8,face="bold"),
           strip.text.x = element_text(size=8, face="bold"),
           legend.text = element_text(size=8),
           legend.title = element_text(size=8))+ guides(fill = guide_legend(override.aes = list(size=5)),
color= "none")


a=d2+facet_grid(~group_status,scale="free",space="free")+theme(strip.text.x = element_blank(),strip.text.y = element_text(
    size =8, face="bold"))


d=ggplot(df2[df2$site == 'Volar forearm', ]) +
  aes(x=subject_id, y=value, fill=Target,color=Target) +
  geom_bar(stat="identity", position="fill", width=0.8) +
  theme_bw(base_size=12) +
  ylab('Relative Abundance/%') +xlab("Status")+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0', '', '50', '', '100'), limits=c(0,1)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=8,face="bold"),
        axis.title.y=element_text(face="bold",size=8), axis.text.y=element_text(face="bold",size=8)) +
  theme(legend.key.size=unit(0.5, "cm"), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(legend.position='right') +
  theme(legend.box="vertical", legend.title=element_text(), legend.text=element_text()) +
  scale_fill_manual("Classification", values=color_pall) +scale_color_manual(values=color_pall)+
  theme(strip.text.x=element_text(), strip.text.y=element_text(angle=0, ), strip.background=element_rect(fill="white")) +
  guides(fill = guide_legend(reverse=FALSE, ncol=1))


d2=d+theme(axis.text.x = element_blank(),
           axis.text.y = element_text(size=8,face="bold"),
           axis.title.x = element_text(size=8,face="bold"),
           axis.title.y = element_text(size=8,face="bold"),
           strip.text.x = element_text(size=8, face="bold"),
           legend.text = element_text(size=8),
           legend.title = element_text(size=8))+ guides(fill = guide_legend(override.aes = list(size=5)),
color= "none")


b=d2+facet_grid(~group_status,scale="free",space="free")+theme(strip.text.x = element_blank(
  ),strip.text.y = element_text(
    size =8, face="bold"))



d=ggplot(df2[df2$site == 'Toeweb', ]) +
  aes(x=subject_id, y=value, fill=Target,color=Target) +
  geom_bar(stat="identity", position="fill", width=0.8) +
  theme_bw(base_size=12) +
  ylab('Relative Abundance/%') +xlab("Status")+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0', '', '50', '', '100'), limits=c(0,1)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_text(size=8,face="bold"),
        axis.title.y=element_text(face="bold",size=8), axis.text.y=element_text(face="bold",size=8)) +
  theme(legend.key.size=unit(0.5, "cm"), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  theme(legend.position='right') +
  theme(legend.box="vertical", legend.title=element_text(), legend.text=element_text()) +
  scale_fill_manual("Classification", values=color_pall) +scale_color_manual(values=color_pall)+
  theme(strip.text.x=element_text(), strip.text.y=element_text(angle=0, ), strip.background=element_rect(fill="white")) +
  guides(fill = guide_legend(reverse=FALSE, ncol=1))


d2=d+theme(axis.text.x = element_blank(),
           axis.text.y = element_text(size=8,face="bold"),
           axis.title.x = element_text(size=8,face="bold"),
           axis.title.y = element_text(size=8,face="bold"),
           strip.text.x = element_text(size=8, face="bold"),
           legend.text = element_text(size=8, face="italic"),
           legend.title = element_text(size=8))+ guides(fill = guide_legend(override.aes = list(size=5)),
color= "none")


c=d2+facet_grid(~group_status,scale="free",space="free")+theme(strip.text.x = element_blank(
  ),strip.text.y = element_text(
    size =8, face="bold"))

q=ggarrange(a+theme(legend.position="none"),b+theme(legend.position="none",
axis.text.y=element_blank(),axis.title.y=element_blank()),c+theme(legend.position="none",
axis.text.y=element_blank(),axis.title.y=element_blank()),
nrow=1,ncol=3,widths=c(1,0.8,0.8))

ggsave("ti_det_compo_2006.svg",q,height=2,width=8)

c=c+theme(legend.position="bottom")

ggsave("ti_det_compo_2006_leg.svg",get_legend(c),height=2,width=8)






