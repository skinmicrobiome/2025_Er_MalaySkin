


library(pacman)

pacman::p_load(readr,plyr, dplyr, ggpubr, vegan, plotly, plyr, ggplot2)


setwd("C:/Users/adikx/OneDrive/Desktop/r_plot/plo1/PUB/new_db")

z1=readr::read_csv("rate_growh.csv")

library(reshape2)

w=melt(z1)


w2=filter(w, variable != "max_rate")

w2$variable=as.numeric(w2$variable)


x=ggplot(data=w2, aes(x=variable, y=value, group=Sample_ID,color=village)) +
  geom_line(linewidth=0.3)+
  theme_classic()+
labs(y="Diameter/cm", x="Time/ weeks",color="Village")+
  theme( axis.text.x = element_text(size=8,angle=90,vjust=0.5,hjust=0.5),
         axis.text.y = element_text(size=8),
         axis.title.x =element_text(size=8),
         axis.title.y = element_text(size=8),legend.position="none",
           plot.title=element_text(size=8))+scale_color_manual(values=c("#004e00","#27ff27","#b1ffb1"))


ggsave("Fig_S1D.svg",x+theme(legend.position="right"),height=4,width=5)


