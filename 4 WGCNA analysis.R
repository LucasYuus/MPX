rm(list = ls())  
options(stringsAsFactors = F)
setwd("")
library(WGCNA)
data <- read.table("", sep = "\t", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
data <- read.csv("", row.names = 1)
load("")
group_list <- as.data.frame(group_list)
data <- dat1
dim(data)
fpkm = data
names(fpkm)
class(data)
data <- as.data.frame(dat)
data$mean = apply(data, 1, mean) # 1 represents each row, 2 represents each column
data = data[order(data$mean, decreasing = TRUE),]
five = data[1:8000,]
five = five[, -492] # This step is necessary to remove the newly created mean column
datExpr0 = t(five)
# Check data
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#### If returns FALSE, correct missing values
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
######### Sample clustering
sampleTree = hclust(dist(datExpr0), method = "average")
plot(sampleTree)
par(mar = c(0, 5, 2, 0)) # Numbers represent the distance from the edge of the plot
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "",
     cex.axis = 0.9, cex.main = 1.2, cex.lab = 0.5, cex = 0.5)
# Optimization
par(cex = 0.1)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
## Remove outlier samples
abline(h = 110, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = 110, minSize = 3)
table(clust)
View(clust)
keepSamples = (clust == 1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) # Generate a new matrix without outlier samples
## Expression data processing complete

##
# Calculate soft threshold
# Load the WGCNA package
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the data saved in the first part

powers = c(c(1:10), seq(from = 12, to = 20, by = 2)) # Cannot be changed
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
# sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 1.5
par(mar = c(5, 5, 3, 1))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
# This line corresponds to using an R^2 cut-off of h
abline(h = 0.85, col = "red") # h is generally greater than 0.85

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

power = sft$powerEstimate # Automatically select the best soft threshold
power
# Construct network
softPower = power
adjacency = adjacency(datExpr0, power = softPower)

# Convert adjacency matrix to TOM matrix
TOM = TOMsimilarity(adjacency);
dissTOM = 1 - TOM
## Cluster topology matrix
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# Adjust clustering distribution
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Match modules with colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
## Merge similar modules
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25###### Cutting height can be modified, 75% similarity
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
table(mergedColors)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## Save results
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "stepBystep.RData")

### Relationship between modules and traits - matching process
allTatits <- read.delim("", sep = "\t", check.names = F, stringsAsFactors = F, header = T)
allTatits <- read.csv("")
dim(allTatits)
names(allTatits)
femalSamples = rownames(datExpr0)
traitRows = match(femalSamples, allTatits$ID)
datTraits = allTatits[traitRows, -1]
rownames(datTraits) = allTatits[traitRows, 1]

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs0[1:6, 1:6]
MEs = orderMEs(MEs0)
MEs[1:6, 1:6]
# Calculate the correlation and significance of each module with each trait
# datTraits = datTraits[,-c(1,2,3,4,5,6,7,8,9,10,11,12,13,17,18,19,20,21,22,23,24,25)] # Customize

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitCor[1:6, 1:6]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples); # Calculate significance
moduleTraitPvalue[1:6, 1:6]

sizeGrWindow(12, 12)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Plot the heatmap
library(lattice)
library(ggplot2)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               main = paste("Module-trait relationships"))

