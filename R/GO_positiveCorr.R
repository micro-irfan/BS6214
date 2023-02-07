library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

high_corr <- "CYP2D6_bestCorrelation.txt"

res_signi <- read.csv(high_corr,header=FALSE)

genes_to_test <- res_signi$V1
genes_to_test
GO_results <- enrichGO(gene = genes_to_test, 
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL", 
                       ont = "BP",
                       pAdjustMethod = "none",
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1,
                       readable = TRUE, 
                       pool = TRUE)

as.data.frame(GO_results)
cluster_summary <- data.frame(GO_results)

write.csv(cluster_summary, "clusterProfiler_CYP2D6.csv")
fit1 <- dotplot(GO_results, showCategory=15)
png("GO_dot_plot_CYP2D6.png", res = 250, width = 1400, height = 1800)
print(fit1)
dev.off()

fit2 <- emapplot(GO_results, showCategory = 50)
png("GO_emap_plot_CYP2D6.png", res = 250, width = 1400, height = 4000)
print(fit2)
dev.off()

fit <- plot(barplot(GO_results, showCategory = 15))
png("GO_plot_CYP2D6.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()
