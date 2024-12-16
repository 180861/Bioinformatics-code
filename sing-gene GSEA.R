library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(data.table)
library(biomaRt)
library(ReactomePA)
library(enrichplot)
genelist_input <- fread(file="GSEA.txt", header = T, sep='\t', data.table = F)
genename <- as.character(genelist_input[,1])
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ]
colnames(gene_map)[1]<-"Gene"
temp<-inner_join(gene_map,genelist_input,by = "Gene")
temp<-temp[,-1]
temp<-na.omit(temp)
temp$logFC<-sort(temp$logFC,decreasing = T)
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=0.05)
kegg_results<-as.data.frame(KEGG_gseresult)
write.csv (kegg_results, file ="KEGG_gseresult.csv")

kk2 <- gseKEGG(geneList= geneList,
               organism     = 'hsa',               
               nPerm        = 10000,               
               minGSSize    = 10,               
               maxGSSize    = 200, 
               pvalueCutoff = 0.05,  
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)
num=5
pdf(paste0("2.","down_GSEA.pdf"),width = 5,height = 5)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$NES),num)])
dev.off()
pdf(paste0("2.","up_GSEA.pdf"),width = 5,height = 5)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$NES),num)])
dev.off()
