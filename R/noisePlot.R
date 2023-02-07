library(tidyr)
library(ggplot2)
library(scales)
library(tidyverse)

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

file       <- "Human_genes.mx"
annotation <- "annotationR.csv"

coldata <- read.csv(annotation,row.names=1,header=TRUE)
data    <- read.csv(file,row.names=1,header=TRUE,sep='\t')
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
data <- data+1

noise <- function(genedata, start, end) {
  noise_list <- c()
  column_list <- c()
  col <- colnames(genedata)
  n <- nrow(genedata)
  for(x in start:(end-1)) {
    tmp = x + 1
    for (y in tmp:end){
      ni <- 2*((genedata[x] - genedata[y])**2)/((genedata[x] + genedata[y])**2)
      
      # Replace NaN values with 0 for list
      ni <- rapply( ni, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
      
      # print (ni[1:50,])
      nav <- (1/n) * (c(sum(ni)))
      column_name <- paste(col[x], '/', col[y], sep="")
      column_list <- append(column_list, column_name)
      noise_list <- append(noise_list, nav)
    }
  }
  combine <- setNames(as.list(noise_list), column_list)
  return (combine)
}

liver <- noise(data, 1, 24)
liver[1]
liver.mean <- mean(unlist(liver))
liver.sd   <- sd(unlist(liver))
liver.mean
liver.sd


kidney <- noise(data, 25, 44)
kidney[1]
kidney.mean <- mean(unlist(kidney))
kidney.sd   <- sd(unlist(kidney))
kidney.mean
kidney.sd

heart <- noise(data, 45, 68)
heart
heart.mean <- mean(unlist(heart))
heart.sd   <- sd(unlist(heart))
heart.mean
heart.sd

adipose <- noise(data, 69, 93)
adipose
adipose.mean <- mean(unlist(adipose))
adipose.sd   <- sd(unlist(adipose))
adipose.mean
adipose.sd

df <- data.frame(Reduce(rbind, adipose))
colnames(df)[1] = "Noise"
df$Organ <- 'Adipose'

df1 <- data.frame(Reduce(rbind, heart))
colnames(df1)[1] = "Noise"
df1$Organ <- 'Heart'

df2 <- data.frame(Reduce(rbind, kidney))
colnames(df2)[1] = "Noise"
df2$Organ <- 'Kidney'

df3 <- data.frame(Reduce(rbind, liver))
colnames(df3)[1] = "Noise"
df3$Organ <- 'Liver'

total <- rbind(df, df1, df2, df3)
p<-ggplot(total, aes(x=Organ, y=Noise))  + geom_jitter(position=position_jitter(0.2), cex=1.2, shape=17)

p + stat_summary(total.y=mean, geom="point", shape=18,
                 size=5, color="red")  

p + stat_summary(fun.data=mean_sdl, mult=1, 
                   geom="pointrange", color="red")



