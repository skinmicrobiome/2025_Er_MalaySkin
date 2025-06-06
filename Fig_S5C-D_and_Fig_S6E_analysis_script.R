
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
library(Maaslin2)


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")


ps3=readRDS("phyloseq_0505.rda")

tax=data.frame(tax_table(ps3))
tax$sp2=paste(tax$Genus,tax$Species)
tax2=tax_table(tax)
taxa_names(tax2)=rownames(tax)

tax_table(ps3)=tax2

colnames(tax_table(ps3))=c("Kingdom","Phylum","Class","Order","Family","Genus",
"Species","sp2")

###Load data
n_taxa = ntaxa(ps3)
n_samp = nsamples(ps3)
# Metadata
meta_data = meta(ps3)
# Taxonomy table
taxonomy = tax_table(ps3)
# Absolute abundances
otu_absolute = abundances(ps3)

#Add taxanomy label to OTU table
taxonomy <- data.frame(tax_table(ps3))
taxonomy <- taxonomy %>% rownames_to_column(var = "taxa_id")
otu_absolute<- data.frame(abundances(ps3))
otu_absolute <- otu_absolute %>% rownames_to_column(var = "taxa_id")
otu_taxa<- dplyr::left_join(otu_absolute, taxonomy)


#set new rowname
rownames(otu_taxa) <- otu_taxa$sp2

#Remove unwanted column (use colnames(otu_taxa) command, include the sample_id column only)
otu_taxa2 <- otu_taxa[,c(2:352)]

#adjust the random effects as u need, or remove it if not necessary),
#adjust output filename if need to


meta_data$Group=as.factor(meta_data$Group)
meta_data$Gender=as.factor(meta_data$Gender)
meta_data$village=as.factor(meta_data$village)
meta_data$subject_id=as.factor(meta_data$subject_id)

fit_data1 = Maaslin2(
  input_data = otu_taxa2, 
  input_metadata = meta_data,
  normalization = "TSS",
  transform ="LOG",
  output = "vill_070335", 
  min_prevalence = 0,
  fixed_effects = c("village"),
  reference=c("village,Kuala Lumpur"),
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  max_significance = 0.05,
  random_effects=c("Age.Group","Gender","site_specific"))



#disease


ps4=subset_samples(ps3, !village %in% c("Washington DC","Kuala Lumpur"))
ps4=subset_samples(ps4, Disease %in% c("Negative","Tinea imbricata"))


sample_data(ps4)$group_status <- ifelse(sample_data(ps4)$Disease == "Negative", "Negative", ifelse(sample_data(ps4)$Disease == "Tinea imbricata", sample_data(ps4)$a_u,NA))


sample_data(ps4)$group_status=recode(sample_data(ps4)$group_status,"A"="With lesion","U"="Without lesion")


n_taxa = ntaxa(ps4)
n_samp = nsamples(ps4)
# Metadata
meta_data = meta(ps4)
# Taxonomy table
taxonomy = tax_table(ps4)
# Absolute abundances
otu_absolute = abundances(ps4)

#Add taxanomy label to OTU table
taxonomy <- data.frame(tax_table(ps4))
taxonomy <- taxonomy %>% rownames_to_column(var = "taxa_id")
otu_absolute<- data.frame(abundances(ps4))
otu_absolute <- otu_absolute %>% rownames_to_column(var = "taxa_id")
otu_taxa<- dplyr::left_join(otu_absolute, taxonomy)


#set new rowname
rownames(otu_taxa) <- otu_taxa$sp2

#Remove unwanted column (use colnames(otu_taxa) command, include the sample_id column only)
otu_taxa2 <- otu_taxa[,c(2:127)]

#adjust the random effects as u need, or remove it if not necessary),
#adjust output filename if need to

meta_data$Group=recode(meta_data$Group,"Teen"="Kids")

meta_data$Category=recode(meta_data$Category,"Obese"="Overweight")


meta_data$Group=as.factor(meta_data$Group)
meta_data$Gender=as.factor(meta_data$Gender)
meta_data$village=as.factor(meta_data$village)
meta_data$subject_id=as.factor(meta_data$subject_id)

fit_data1 = Maaslin2(
  input_data = otu_taxa2, 
  input_metadata = meta_data,
  normalization = "TSS",
  transform ="LOG",
  output = "ti_det_160324_kids", 
  min_prevalence = 0.2,
  fixed_effects = c("group_status"),
  reference=c("group_status,Negative"),
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  max_significance = 0.05,
  random_effects=c("Age.Group","Gender","village","site_specific"))





