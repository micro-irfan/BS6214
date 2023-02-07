library(tidyr)
library(ggplot2)
library(scales)
library(tidyverse)
library(reshape2)
library("data.table")

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

file       <- "Human_genes.tpm.csv"
annotation <- "annotationR.csv"

coldata <- read.csv(annotation, row.names=1,header=TRUE)
coldata$run <-  rownames(coldata)

data <- read.csv(file,row.names=1,header=TRUE,sep=',')
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]

## Distribution plot for all sample

data$refseq <- rownames(data)

# coldata[coldata$run == 'SRR2087367.out.bam',]$expriment

data.plot <- melt(data)   
data.plot$tissue <- coldata$expriment[match(data.plot$variable,coldata$run)]

p <- ggplot(data.plot)
p <- p + theme(legend.position = "none") + xlab('TPM') + ylab('DENSITY')
p <- p + geom_density(aes(x=value, colour=variable)) +
         scale_colour_manual(values = coldata$colour) +
         scale_x_continuous(trans = 'log10',
                            limits = c(10^-3, 10^4),
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))  

file_name <- paste0("TPM-Density-AllCodingGenes", ".pdf")
ggsave(p, file = file_name)

p1 <- ggplot(data.plot)
p1 <- p1 + theme(legend.position = "none") + xlab('TPM') + ylab('DENSITY')
p1 <- p1 + geom_density(aes(x=value, colour=tissue)) +
         scale_colour_manual(values = c('green', 'purple', 'blue', 'red')) +
         scale_x_continuous(trans = 'log10',
                            limits = c(10^-3, 10^4),
                            breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))  

file_name <- paste0("TPM-Density-AllCodingGenes.grouped", ".pdf")
ggsave(p1, file = file_name)

file       <- "Human_genes.pharmacogenes.tpm.csv"
annotation <- "annotationR.csv"

coldata <- read.csv(annotation, row.names=1,header=TRUE)
coldata$run <-  rownames(coldata)

data <- read.csv(file,row.names=1,header=TRUE,sep=',')
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]

## Distribution plot for all sample

data$refseq <- rownames(data)

# coldata[coldata$run == 'SRR2087367.out.bam',]$expriment

data.plot <- melt(data)   
data.plot$tissue <- coldata$expriment[match(data.plot$variable,coldata$run)]

p <- ggplot(data.plot)
p <- p + theme(legend.position = "none") + xlab('TPM') + ylab('DENSITY')
p <- p + geom_density(aes(x=value, colour=variable)) +
  scale_colour_manual(values = coldata$colour) +
  scale_x_continuous(trans = 'log10',
                     limits = c(10^-3, 10^4),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))  

file_name <- paste0("TPM-Density-Pharmacogenes", ".pdf")
ggsave(p, file = file_name)

p1 <- ggplot(data.plot)
p1 <- p1 + theme(legend.position = "none") + xlab('TPM') + ylab('DENSITY')
p1 <- p1 + geom_density(aes(x=value, colour=tissue)) +
  scale_colour_manual(values = c('green', 'purple', 'blue', 'red')) +
  scale_x_continuous(trans = 'log10',
                     limits = c(10^-3, 10^4),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))  

file_name <- paste0("TPM-Density-Pharmacogenes.grouped", ".pdf")
ggsave(p1, file = file_name)

