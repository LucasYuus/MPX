rm(list = ls())
library(tidyverse)
library(readxl)
load(file = "")
exp <- ssgsea
exp <- as.data.frame(exp)
setwd("")  

genes <- rownames(exp)

exp <- t(exp)  
exp <- data.table::fread("", data.table = F)
exp <- read.csv("", row.names = 1)

rownames(exp) = exp[,1]
exp = exp[,-1]
exp <- exp %>% column_to_rownames("") %>% scale() %>% t() %>% as.data.frame()

rt <- read_excel(path = "")

genes <- rt$`cholesterol metabolism related gene`

hc_geneExp <- exp[genes, ]
hc_geneExp <- as.matrix(hc_geneExp)
hc_geneExp <- exp

library(ConsensusClusterPlus)
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)

Cluster <- ConsensusClusterPlus(d = hc_geneExp,  # Analysis matrix
                                maxK = 9,  # Maximum number of clusters
                                reps = 1000,  # Resampling times
                                pItem = 0.8,  # Resampling proportion of samples
                                clusterAlg = "km",  # Clustering algorithm used; options include "hc" (hclust), "pam", "km" (k-means)
                                innerLinkage = "ward.D2", 
                                finalLinkage = "ward.D2",
                                distance = "euclidean",  
                                seed = 44,  # Set random seed
                                title = "Consensus Cluster") 

#
Kvec <- 2:4
x1 = 0.1; x2 = 0.9        
PAC <- rep(NA, length(Kvec)) 
names(PAC) <- paste("K=", Kvec, sep="")  
for(i in Kvec){
  M = Cluster[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])          # Calculate the consensus matrix
  PAC[i-1] = Fn(x2) - Fn(x1)
} 

optK <- Kvec[which.min(PAC)]  # Optimal K value
optK
#

annCol <- data.frame(Cluster = paste0("Cluster",
                                      Cluster[[2]][["consensusClass"]]),
                     row.names = colnames(hc_geneExp))
table(annCol$Cluster)
write.csv(annCol, "")

identical(rownames(annCol), colnames(hc_geneExp))

# PCA dimensionality reduction
#BiocManager::install("factoextra")
#install.packages("estimability")
library(FactoMineR)
library(factoextra)

exp <- as.data.frame(t(hc_geneExp))
dat <- PCA(exp, graph = FALSE)
pca <- fviz_pca_ind(dat,
                    col.var = "steelblue",
                    title = "PCA - Biplot",
                    geom.ind = "point",
                    col.ind = annCol$Cluster,
                    addEllipses = TRUE,
                    pointsize = 2.5,  # Size of points
                    ellipse.type = "convex",  # Type of outer boundary, convex polygon
                    axes.linetype = "blank",  # Hide zero lines
                    legend.title = "Groups")
pca
dev.off()