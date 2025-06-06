setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db/input_file_r_code_coded")

library(pacman)

pacman::p_load(dplyr,phyloseq,tidyverse, microbiome,treeio,ggtree,ape)



tree2=read.tree("cory_recoded_tree.nwk")
df2=read.csv("tree_name.csv")

ordering_index <- match(df2$recoded, unique(tree2$tip.label))

# Use the ordering index to reorder df2
df3 <- df2[order(ordering_index),]

tree2$tip.label=df3$Species

circ=ggtree(tree2) +
   theme_tree2() + 
  geom_tiplab(align=TRUE, linetype='dashed', as_ylab=TRUE,linesize=.3,color='black',size=8)


ggsave("Fig_5A",circ,height=6,width=7.7,dpi=320)
