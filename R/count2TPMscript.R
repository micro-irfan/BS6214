library("DESeq2")

setwd("C:/Irfan/Personal/Masters/BS6214/Assignment1/src")

file        <- "Human_genes.mx"
file_length <- "Human_genes_length.mx"

data        <- read.csv(file,row.names=1,header=TRUE,sep='\t')
data_length <- read.csv(file_length,row.names=1,header=TRUE,sep='\t')

Counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

tpm <- Counts_to_tpm(data, data_length$Length)
tpm
write.csv(tpm, "Human_genes.tpm.csv", row.names=TRUE)
