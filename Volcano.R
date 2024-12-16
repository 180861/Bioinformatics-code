library(ggpubr)
library(ggthemes)
df <- read.csv("OB-Vol.csv",header=T,stringsAsFactors = F)
df$logP <- -log10(df$pvalue)
df$group[((df$log2FoldChange<1)&(df$log2FoldChange>=-1))|(df$pvalue>=0.05)] <- "not significant"
df$group[(df$log2FoldChange>=1)&(df$pvalue<0.05)] <- "up-regulated"
df$group[(df$log2FoldChange<=-1)&(df$pvalue<0.05)] <- "down-regulated"
table(df$group)

ggscatter(df,x="log2FoldChange",y="logP",color = "group",
          palette=c("#00ba38","gray","#f8766d"),size=1,
          font.label = 8,xlab = "log2FoldChange",
          ylab="-log10(P value)")+theme_base()+
  xlim(-4,4)+ylim(0,15)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = c(-1,1),linetype="dashed")
