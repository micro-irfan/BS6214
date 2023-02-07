

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

file <- "Human_genes.liver.mx"
data <- read.csv(file,row.names=1,header=TRUE,sep='\t')
data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]

data.pca <- prcomp(data,
               center = TRUE,
               scale. = TRUE)
summary(data.pca)
data.pca

library(ggfortify)
data.pca.plot <- autoplot(data.pca,
                          data = data,)
data.pca.plot

# biplot.data.pca <- biplot(data.pca)
# biplot.data.pca

library(factoextra)

res.pca <- prcomp(data, scale = TRUE)
fviz_eig(res.pca)

ind <- get_pca_ind(res.pca)
test <- as.data.frame(ind$contrib)

test <- test[order(test$Dim.1,decreasing=TRUE),]
top100 <- rownames(test)[1:100]
lapply(top100, write, "LiverTop100PCA.txt", append=TRUE, ncolumns=1000)


test[c('ENSG00000140505.6'),]

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
