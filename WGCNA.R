library(WGCNA)
options(stringsAsFactors = FALSE)
femData = read.csv("Merge-OB.csv") 
dim(femData)
names(femData) 
datExpr0 = as.data.frame(t(femData[, -c(1)])) 
names(datExpr0) = femData$geneNames
rownames(datExpr0) = names(femData)[-c(1)]
View(head(datExpr0))

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


vars_res <- apply(datExpr0, 2, var)
per_res <- quantile(vars_res, probs = seq(0, 1, 0.25))
per_res
upperGene <- datExpr0[, which(vars_res > per_res[4])]
dim(upperGene)


sampleTree = hclust(dist(upperGene), method = "average")
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


abline(h = 60, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)   

keepSamples = (clust==1)  
datExpr = upperGene[keepSamples, ] 
nGenes = ncol(datExpr)
nGenes
nSamples = nrow(datExpr)
nSamples


traitData = read.csv("ClinicalTraits.csv")
femaleSamples = rownames(datExpr) 
traitRows = match(femaleSamples, traitData$ID)
datTraits = traitData[traitRows, -1] 
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()


sampleTree2 = hclust(dist(datExpr), method = "average") 
traitColors = numbers2colors(datTraits, signed = FALSE) 
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")



enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2)) 
cex1 = 0.85
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")  

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


sft$powerEstimate

net <- blockwiseModules(datExpr, power = 5, 
                        maxBlockSize = nGenes, TOMType = "unsigned", 
                        minModuleSize = 30, reassignThreshold = 0, 
                        mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = F, verbose = 3)


table(net$colors)



sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]


nGenes = ncol(datExpr)

nSamples = nrow(datExpr) 
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) 
moduleTraitCor = cor(MEs, datTraits, use = "p") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) 
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 7, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



OB = as.data.frame(datTraits$OB)
names(OB) = "OB"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, OB, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(OB), sep="")
names(GSPvalue) = paste("p.GS.", names(OB), sep="")


probes = names(datExpr) 
geneInfo0 = data.frame(probes = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue);

modOrder = order(-abs(cor(MEs, OB, use = "p")))
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor); 
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")
