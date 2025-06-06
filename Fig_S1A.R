setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")
library(pacman)

pacman::p_load(dplyr,phyloseq,tidyverse, microbiome,treeio,ggtree,ape)

tree2=read.tree("tricho.iqtree_0305_coded.nwk")
df2=read.csv("tri_tree.csv")


# Create an ordering index based on the match between df2$id and unique(tree2$tip.label)
ordering_index <- match(df2$id2, unique(tree2$tip.label))

# Use the ordering index to reorder df2
df3 <- df2[order(ordering_index),]


tree2$tip.label=df3$tip

df_v <- data.frame(Location=df3$Location
)



rownames(df_v) <- tree2$tip.label

df_ph<- data.frame(Phylum=df3$Host
)

rownames(df_ph) <- tree2$tip.label




circ=ggtree(tree2) +
   theme_tree2() + 
  geom_tiplab(align=TRUE, linetype='dashed', as_ylab=TRUE,linesize=.3,color='black',size=8,face="italic")

library(ggsci)


p1 <- gheatmap(circ, df_ph, offset=.01, width=.2,
               colnames_angle=45, colnames_offset_y = -1.1,hjust=1) +
  scale_fill_manual(values=c("lightsalmon2","seagreen"))+
  labs(fill="Host")
colors = c("#FFBC00","#FF8000","#C6DBAF",
                          "#6BBF59","#21D375","#08A045","#25A18E")

library(ggnewscale)
p2 <- p1 + new_scale_fill()
p3=gheatmap(p2, df_v, offset=.05, width=.2,
            colnames_angle=45, colnames_offset_y = -1.1,hjust=1) +
  scale_fill_manual(values=c("#27ff27","#25A18E","#b1ffb1","honeydew","#004e00"))+labs(fill="Source")+
  theme( legend.text=element_text(size=8,face="bold"),
         legend.title=element_text(size=8,face="bold")
  )



ggsave("tri_0103.svg",p3+theme(legend.position="none"),height=5,width=4)


