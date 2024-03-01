

library(ggplot2)

my_theme<- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.border = element_rect(colour="black",fill=NA))+
  theme(axis.text.x = element_text(size=8, color = "black"),  
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.x = element_text(size=8, color = "black"),
        axis.title.y= element_text(size=8,  color = "black"),
        legend.title = element_text(size=8,  color = "black"))+
  theme(panel.grid =element_blank())  +
  theme(axis.line = element_line(size=0.3, colour = "black")) +
  theme(panel.grid.major=element_line(colour=NA), panel.border = element_blank())+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"),angle=0),
        # axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        legend.background = element_rect(fill="white"),
        legend.key.size = unit(15, "pt"),
        legend.key.height = unit(15, "pt"),
        legend.key.width = unit(15, "pt"))


bar_plot<-ggplot(data=data, mapping=aes(x=group,y=pCR))+
  geom_bar(stat="identity",fill='#4091c9')+
  my_theme+
  labs(y="pCR rates",x="")+
  geom_text(aes(x=group,
                y = pCR+7, 
                label = round(pCR,1)),   
            color="black") +
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"),angle=60))
bar_plot



