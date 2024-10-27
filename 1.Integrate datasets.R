
library(limma)
library(sva)
library(tinyarray)
setwd("")
# File names
files=c("")      

geneList=list()
for(i in 1:length(files)){
	fileName=files[i]
    rt=read.table(fileName, header=T, sep="\t", check.names=F)
    header=unlist(strsplit(fileName, "\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
	fileName=files[i]
    header=unlist(strsplit(fileName, "\\.|\\-"))
    rt=read.table(fileName, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    
}

# Integrate data
sample=read.table("", header=F, sep="\t", check.names=F, row.names = 1)
allTab=allTab[,rownames(sample)]
allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file="", sep="\t", quote=F, col.names=F)

# batchType: Source of series; modeType: Sample grouping type
batchType=as.character(as.data.frame(strsplit(rownames(sample), "_"))[1,])
modeType=as.vector(sample[,1])
mod=model.matrix(~as.factor(modeType))
normalizeTab=ComBat(allTab, batchType, mod, par.prior=TRUE)
normalizeTabout=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTabout, file="", sep="\t", quote=F, col.names=F)

############# Before correction
# Sample grouping
type=factor(x = colnames(allTab), labels =modeType)
summary(type)
pdf(file = "", width=5.5, height = 4.5)
draw_pca(allTab, type)
dev.off()
# Source of series for samples
type=factor(x = colnames(allTab), labels =batchType)
summary(type)
pdf(file = "", width=5.5, height = 4.5)
draw_pca(allTab, type)
dev.off()
############# After correction
# Sample grouping
type=factor(x = colnames(normalizeTab), labels =modeType)
summary(type)
pdf(file = "", width=5.5, height = 4.5)
draw_pca(normalizeTab, type)
dev.off()
# Source of series for samples
type=factor(x = colnames(normalizeTab), labels =batchType)
summary(type)
pdf(file = "", width=5.5, height = 4.5)
draw_pca(normalizeTab, type)
dev.off()