library(tidyverse)
library(data.table)
library(GSVA)
library(ggsci)
library(tidyr)
library(ggpubr)

cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"
type <- split(cellMarker,cellMarker$celltype)
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})
expr <- data.table::fread("exp2.txt",data.table = F)   
rownames(expr) <- expr[,1]   
expr <- expr[,-1]
expr <- log(expr+1) 
expr <- as.matrix(expr)

gsva_data <- gsva(expr,cellMarker, method = "ssgsea")


data <- gsva_data %>% t() %>% as.data.frame()
group_list = ifelse((str_sub(row.names(data),14,15))<10,'1','0') %>% as.data.frame()
colnames(group_list)[1] <- 'group'
data$group <- group_list$group
data <- data %>% rownames_to_column("sample")
dat <- gather(data,key=ssGSEA,value = Expression,-c(group,sample))

ggboxplot(dat, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1))

