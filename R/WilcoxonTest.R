suppressWarnings(library(edgeR, quietly = T))

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

file       <- "Human_genes.mx"
annotation <- "annotationR.csv"

conditions <- read.table(annotation, row.names=1, header=TRUE, sep=',')
readCount  <- read.table(file=file, header = T, row.names = 1, stringsAsFactors = F,check.names = F)

liverConditions<-replace(conditions$expriment, conditions$expriment!='Liver', 'Not Liver')

y <- DGEList(counts=readCount,group=liverConditions)
##Remove rows consistently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~liverConditions, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")

conditionsLevel<-levels(liverConditions)
dataCon1=count_norm[,c(which(liverConditions=='Not Liver'))]
dataCon2=count_norm[,c(which(liverConditions=='Liver'))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
foldChanges

outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)

summarize<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
nrow(summarize)
upLFC<-subset(summarize, (log2foldChange > 2 & pValues < 0.05 & FDR<0.05))
nrow(upLFC)
downLFC<-subset(summarize, (log2foldChange < -2 & pValues < 0.05 & FDR<0.05))
nrow(downLFC)

rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05
write.table(outRst[outRst$FDR<fdrThres,], file="Liver.WilcoxonTest.rst.tsv",sep="\t", quote=F,row.names = T,col.names = T)

outRst

library(EnhancedVolcano)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- data.frame(log2FoldChange=foldChanges, pValues=pvalues, FDR=fdr)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(pValues), xlim = c(-12, 12), ylim = c(0,20), pch=20, main="Liver Transcriptomic (Wilcoxon)", cex=0.75, xlab=bquote(~Log[2]~"(fold-change)"), ylab=bquote(~-log[10]~"(P-value)")))

with(subset(topT, pValues<0.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pValues), pch=20, col="red", cex=0.5))

#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05

abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-0.75, col="black", lty=3, lwd=2.0)
abline(v=.75, col="black", lty=3, lwd=2.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pValues[topT$pValues<0.01], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

