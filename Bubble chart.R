library(ggplot2)
library(Cairo)
data <- read.table("Pathway.txt",header = T,sep = "\t")
result <- ggplot(data,aes(RichFactor,reorder(Pathway,1/FDR)))+
  geom_point(aes(size=as.factor(GeneNum),color=FDR))+
  scale_colour_gradient(low="red",high="blue")+
  labs(color=expression(FDR),size="Gene Number",x="Rich Factor",y="Pathway")+
  theme_bw()+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(aspect.ratio=2,axis.text.x=element_text(size=8,face = "bold",color = "black"),axis.text.y = element_text(size=9,face = "bold",color = "black"),axis.title.x = element_text(size=12,face = "bold",color = "black"),axis.title.y = element_text(size=12,face = "bold",color = "black"),legend.title = element_text(size=10,face = "bold"),legend.text = element_text(size=8,face = "bold"))
result