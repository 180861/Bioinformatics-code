library(GEOquery)
gset <- getGEO("GSE151839",destdir = ".",AnnotGPL = F,getGPL = F)
gset[["GSE151839_series_matrix.txt.gz"]]@annotation

gpl <- data.table::fread("GPL570-55999.txt",sep="\t",header = TRUE)
GPL <- gpl[gpl$'Gene Symbol'!='',]


exprSet <- exprs(gset[[1]])
exprSet <- as.data.frame(exprSet)


a3=data.frame(GPL$ID,GPL$'Gene Symbol')
colnames(a3) <- c("ID","gene.all")

a3$gene.all <- sub("///.*","",a3$gene.all)
a3$gene.all <- sub("- *","",a3$gene.all)

ids=merge(exprSet,a3,by.x=0,by.y="ID")
dim(ids)


ids=ids[!duplicated(ids$gene.all),]
row.names(ids)=ids$gene.all 
ids=ids[,-1]
ids=ids[,-ncol(ids)] 

range(ids)
boxplot(ids[,1:20])

library(limma)
exp=log2(ids+1)
range(exp)

boxplot(exp[,1:20])
exp=normalizeBetweenArrays(exp)
write.csv(exp,"GSE151839.csv",quote = F)
