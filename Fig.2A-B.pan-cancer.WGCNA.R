library(WGCNA)
library(limma)

#------------------------------------ACC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\ACC")

mydata <- read.table("lncRNA.ACC.txt",header=T,sep="\t",check.names=F,row.names=1)

group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata)
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=F,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.ACC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.ACC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"ACC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.ACC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------BLCA----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\BLCA")

mydata <- read.table("lncRNA.BLCA.txt",header=T,sep="\t",check.names=F,row.names=1)

group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata)
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=F,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.BLCA.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.BLCA.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"BLCA.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.BLCA.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------BRCA----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\BRCA")

mydata <- read.table("lncRNA.BRCA.txt",header=T,sep="\t",check.names=F,row.names=1)

group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.BRCA.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.BRCA.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"BRCA.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.BRCA.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------CESC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\CESC")

mydata <- read.table("lncRNA.CESC.txt",header=T,sep="\t",check.names=F,row.names=1)

group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.CESC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.CESC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"CESC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.CESC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------CHOL----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\CHOL")

mydata <- read.table("lncRNA.CHOL.txt",header=T,sep="\t",check.names=F,row.names=1)

group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.CHOL.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.CHOL.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"CHOL.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.CHOL.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------COAD----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\COAD")

mydata <- read.table("lncRNA.COAD.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.COAD.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.COAD.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"COAD.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.COAD.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




# #------------------------------------DLBC----------------------------------
# setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\DLBC")
# 
# mydata <- read.table("lncRNA.DLBC.txt",header=T,sep="\t",check.names=F,row.names=1)
# group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# mydata=mydata[,group==0]
# mydata=t(mydata) 
# #rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
# mydata=t(avereps(mydata))
# #write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)
# 
# datExpr0 = data.frame(t(mydata))
# colnames(datExpr0) <- rownames(mydata)
# rownames(datExpr0) <- colnames(mydata)
# datExpr1<-datExpr0
# m.vars=apply(datExpr0,2,var)
# expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
# datExpr1<-data.matrix(expro.upper)
# gsg = goodSamplesGenes(datExpr1, verbose = 3);
# gsg$allOK
# if (!gsg$allOK){
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0) 
#     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0) 
#     printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
# }
# sampleTree = hclust(dist(datExpr1), method = "average")
# par(cex = 0.7);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# datExpr = as.data.frame(datExpr1)
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 
# pdf("1Threshold.DLBC.pdf",width = 10, height = 5)
# par(mfrow = c(1,2))
# cex1 = 0.9
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence")) +
#   text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#        labels=powers,cex=cex1,col="red")+
#   abline(h=0.85,col="red")
# 
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity")) +
#   text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()
# net = blockwiseModules(datExpr, power = 2,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        #saveTOMs = TRUE,
#                        #saveTOMFileBase = "MyTOM",
#                        verbose = 3)
# 
# table(net$colors)
# mergedColors = labels2colors(net$colors)
# 
# pdf("2module.DLBC.pdf",width = 10, height = 5)
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# MEs = net$MEs;
# geneTree = net$dendrograms[[1]]
# text <- unique(moduleColors)
# for (i  in 1:length(text)) {
#   y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
#   write.csv(y,paste(text[i],"DLBC.csv",sep = "."),quote = F)
# }
# 
# samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)
# 
# sameSample=intersect(row.names(samples), row.names(datExpr))
# samples=samples[sameSample,]
# datExpr=datExpr[sameSample,]
# 
# 
# moduleLabelsAutomatic = net$colors
# moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# moduleColorsWW = moduleColorsAutomatic
# MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
# MEsWW = orderMEs(MEs0)
# modTraitCor = cor(MEsWW, samples, use = "p")
# colnames(MEsWW)
# modlues=MEsWW
# 
# modTraitP = corPvalueStudent(modTraitCor, nSamples)
# textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
# dim(textMatrix) = dim(modTraitCor)
# 
# pdf("3Module-trait.DLBC.pdf",width = 6, height = 6)
# labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
#                xColorWidth = 0.03,
#                ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
#                textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
#                , main = paste("Module-trait relationships"))
# dev.off()




#------------------------------------ESCA----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\ESCA")

mydata <- read.table("lncRNA.ESCA.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.ESCA.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.ESCA.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"ESCA.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.ESCA.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------GBM----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\GBM")

mydata <- read.table("lncRNA.GBM.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.GBM.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.GBM.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"GBM.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.GBM.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------HNSC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\HNSC")

mydata <- read.table("lncRNA.HNSC.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.HNSC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.HNSC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"HNSC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.HNSC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------KICH----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\KICH")

mydata <- read.table("lncRNA.KICH.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.KICH.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.KICH.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"KICH.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.KICH.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------KIRC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\KIRC")

mydata <- read.table("lncRNA.KIRC.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.KIRC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.KIRC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"KIRC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.KIRC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------KIRP----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\KIRP")

mydata <- read.table("lncRNA.KIRP.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.KIRP.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.KIRP.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"KIRP.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.KIRP.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




# #------------------------------------LAML----------------------------------
# setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\LAML")
# 
# mydata <- read.table("lncRNA.LAML.txt",header=T,sep="\t",check.names=F,row.names=1)
# group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
# group=sapply(strsplit(group,""), "[", 1)
# group=gsub("2", "1", group)
# mydata=mydata[,group==0]
# mydata=t(mydata) 
# #rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
# mydata=t(avereps(mydata))
# #write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)
# 
# datExpr0 = data.frame(t(mydata))
# colnames(datExpr0) <- rownames(mydata)
# rownames(datExpr0) <- colnames(mydata)
# datExpr1<-datExpr0
# m.vars=apply(datExpr0,2,var)
# expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
# datExpr1<-data.matrix(expro.upper)
# gsg = goodSamplesGenes(datExpr1, verbose = 3);
# gsg$allOK
# if (!gsg$allOK){
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0) 
#     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0) 
#     printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
# }
# sampleTree = hclust(dist(datExpr1), method = "average")
# par(cex = 0.7);
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# datExpr = as.data.frame(datExpr1)
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 
# pdf("1Threshold.LAML.pdf",width = 10, height = 5)
# par(mfrow = c(1,2))
# cex1 = 0.9
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence")) +
#   text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#        labels=powers,cex=cex1,col="red")+
#   abline(h=0.85,col="red")
# 
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity")) +
#   text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()
# net = blockwiseModules(datExpr, power = 4,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        #saveTOMs = TRUE,
#                        #saveTOMFileBase = "MyTOM",
#                        verbose = 3)
# 
# table(net$colors)
# mergedColors = labels2colors(net$colors)
# 
# pdf("2module.LAML.pdf",width = 10, height = 5)
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# MEs = net$MEs;
# geneTree = net$dendrograms[[1]]
# text <- unique(moduleColors)
# for (i  in 1:length(text)) {
#   y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
#   write.csv(y,paste(text[i],"LAML.csv",sep = "."),quote = F)
# }
# 
# samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)
# 
# sameSample=intersect(row.names(samples), row.names(datExpr))
# samples=samples[sameSample,]
# datExpr=datExpr[sameSample,]
# 
# 
# moduleLabelsAutomatic = net$colors
# moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# moduleColorsWW = moduleColorsAutomatic
# MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
# MEsWW = orderMEs(MEs0)
# modTraitCor = cor(MEsWW, samples, use = "p")
# colnames(MEsWW)
# modlues=MEsWW
# 
# modTraitP = corPvalueStudent(modTraitCor, nSamples)
# textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
# dim(textMatrix) = dim(modTraitCor)
# 
# pdf("3Module-trait.LAML.pdf",width = 6, height = 6)
# labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
#                xColorWidth = 0.03,
#                ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
#                textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
#                , main = paste("Module-trait relationships"))
# dev.off()
# 
# 


#------------------------------------LGG----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\LGG")

mydata <- read.table("lncRNA.LGG.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.LGG.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.LGG.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"LGG.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.LGG.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------LIHC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\LIHC")

mydata <- read.table("lncRNA.LIHC.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.LIHC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.LIHC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"LIHC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.LIHC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------LUAD----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\LUAD")

mydata <- read.table("lncRNA.LUAD.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.LUAD.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.LUAD.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"LUAD.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.LUAD.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------LUSC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\LUSC")

mydata <- read.table("lncRNA.LUSC.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.LUSC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.LUSC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"LUSC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.LUSC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------MESO----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\MESO")

mydata <- read.table("lncRNA.MESO.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.MESO.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.MESO.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"MESO.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.MESO.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------OV----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\OV")

mydata <- read.table("lncRNA.OV.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.OV.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.OV.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"OV.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.OV.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------PAAD----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\PAAD")

mydata <- read.table("lncRNA.PAAD.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.PAAD.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.PAAD.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"PAAD.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.PAAD.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------PCPG----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\PCPG")

mydata <- read.table("lncRNA.PCPG.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.PCPG.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.PCPG.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"PCPG.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.PCPG.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------PRAD----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\PRAD")

mydata <- read.table("lncRNA.PRAD.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.PRAD.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.PRAD.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"PRAD.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.PRAD.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------READ----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\READ")

mydata <- read.table("lncRNA.READ.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.READ.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.READ.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"READ.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.READ.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------SARC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\SARC")

mydata <- read.table("lncRNA.SARC.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.SARC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.SARC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"SARC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.SARC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------SKCM----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\SKCM")

mydata <- read.table("lncRNA.SKCM.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.SKCM.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.SKCM.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"SKCM.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.SKCM.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------STAD----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\STAD")

mydata <- read.table("lncRNA.STAD.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.STAD.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.STAD.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"STAD.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.STAD.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------TGCT----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\TGCT")

mydata <- read.table("lncRNA.TGCT.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.TGCT.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.TGCT.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"TGCT.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.TGCT.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------THCA----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\THCA")

mydata <- read.table("lncRNA.THCA.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.THCA.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.THCA.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"THCA.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.THCA.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------THYM----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\THYM")

mydata <- read.table("lncRNA.THYM.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.THYM.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.THYM.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"THYM.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.THYM.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------UCEC----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\UCEC")

mydata <- read.table("lncRNA.UCEC.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.UCEC.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.UCEC.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"UCEC.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.UCEC.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------UCS----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\UCS")

mydata <- read.table("lncRNA.UCS.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.UCS.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.UCS.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"UCS.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.UCS.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




#------------------------------------UVM----------------------------------
setwd("D:\\BaiduNetdiskWorkspace\\05CAF\\08WGCNA\\UVM")

mydata <- read.table("lncRNA.UVM.txt",header=T,sep="\t",check.names=F,row.names=1)
group=sapply(strsplit(colnames(mydata),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
mydata=mydata[,group==0]
mydata=t(mydata) 
#rownames(mydata)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(mydata))  #用替换的方法删除-01.02,03等等
mydata=t(avereps(mydata))
#write.table(mydata,file="data01.txt",sep="\t",row.names=T,quote=F)

datExpr0 = data.frame(t(mydata))
colnames(datExpr0) <- rownames(mydata)
rownames(datExpr0) <- colnames(mydata)
datExpr1<-datExpr0
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[3])]  #4ǰ25%??3ǰ50%??2ǰ75%??1=ȫ??
datExpr1<-data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("1Threshold.UVM.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)

table(net$colors)
mergedColors = labels2colors(net$colors)

pdf("2module.UVM.pdf",width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"UVM.csv",sep = "."),quote = F)
}

samples=read.csv('normalize-lie.z-score.txt',sep = '\t',row.names = 1)

sameSample=intersect(row.names(samples), row.names(datExpr))
samples=samples[sameSample,]
datExpr=datExpr[sameSample,]


moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW

modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf("3Module-trait.UVM.pdf",width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01,
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()




