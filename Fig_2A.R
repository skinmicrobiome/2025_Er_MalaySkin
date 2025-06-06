
setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

library(pacman)


pacman::p_load(dplyr,phyloseq,tidyverse, microbiome,treeio,ggtree,ape,ggnewscale, cowplot)




tree2=read.tree("recoded_tree_0903.nwk")

df2=read.csv("tree_1602_3.csv")



tree2=keep.tip(tree2, df2$tip2)

ordering_index <- match(df2$tip2, unique(tree2$tip.label))

# Use the ordering index to reorder df2
df3 <- df2[order(ordering_index),]


df_v <- data.frame(village=df3$village
)


rownames(df_v) <- tree2$tip.label


df_x <- data.frame(family=df3$Family
)




rownames(df_x) <- tree2$tip.label






#circ <- ggtree(tree2, layout ="circular")
circ=ggtree(tree2) +
   theme_tree2() + 
  geom_tiplab(align=TRUE, linetype='dashed', as_ylab=TRUE,linesize=.3,color='black',size=8,face="bold")

library(ggsci)

p1 <- gheatmap(circ, df_v, offset=.08, width=.2,
               colnames_angle=45, colnames_offset_y = -1.1,hjust=1) +
  scale_fill_manual(values=c("#25A18E","#004e00","#27ff27","#b1ffb1"))+
  labs(fill="Sources")+ geom_treescale(fontsize = 3)


p4 <- p1 + new_scale_fill()
p5=gheatmap(p4, df_x, offset=0.15, width=.2,
            colnames_angle=45, colnames_offset_y = -1.1,hjust=1) +
  scale_fill_manual(values=c("red","blue","purple","honeydew"))+labs(fill="Family")+
  theme( legend.text=element_text(size=8,face="bold"),
         legend.title=element_text(size=8,face="bold")
  )+ 
  guides(colour = guide_legend(override.aes = list(size=20)))



ggsave("Fig_2A.pdf",p5,height=6,width=7.7,dpi=320)

ggsave("Fig_2A_leg.pdf",get_legend(p5),height=6,width=4,dpi=320)

