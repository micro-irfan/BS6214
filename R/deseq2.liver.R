library("DESeq2")
library(dplyr)

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

file       <- "Human_genes.mx"
annotation <- "annotationR.csv"

conditions <- read.csv(annotation,row.names=1,header=TRUE)
liverConditions<-replace(conditions$expriment, conditions$expriment!='Liver', 'NotLiver')

liverConditions <- as.data.frame(liverConditions)
row.names(liverConditions) <- c(rownames(conditions))

data <- read.csv(file,row.names=1,header=TRUE,sep='\t')
data <- data %>% 
  mutate(Total_ag = rowSums(across(where(is.numeric), 
                                   ~ replace(.x, .x < 10 , NA)), na.rm = TRUE))

data<-subset(data, Total_ag > 10)
drops <- c("Total_ag")
data <- data[ , !(names(data) %in% drops)]

#data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
#data <- data+1

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData   = liverConditions,
                              design    = ~liverConditions)

dds <- DESeq(dds)
res <- results(dds, contrast=c("liverConditions","Liver","NotLiver"), name="test", alpha = 0.05)
res <- lfcShrink(dds, contrast=c("liverConditions","Liver","NotLiver"), res=res, type='ashr')  
write.table(res, file="Liver.DESeq2.csv")

plotMA(res,ylim=c(-8,8))
mcols(res)
summary(res)

library(EnhancedVolcano)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

topT <- as.data.frame(res)

#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), xlim = c(-12, 12), ylim = c(0,305), pch=20, main="Liver Transcriptomic", cex=0.75, xlab=bquote(~Log[2]~"(fold-change)"), ylab=bquote(~-log[10]~"(P-value)")))

with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05

abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-0.75, col="black", lty=3, lwd=2.0)
abline(v=.75, col="black", lty=3, lwd=2.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.01], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
