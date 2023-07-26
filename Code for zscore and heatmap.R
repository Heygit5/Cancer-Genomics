

data<- read.table(file="C:/Users/Anuja/Documents/cancergenomics.txt", 
                  header=TRUE, row.names=1, sep = '\t')
data= data[,-1]
head(data)


sizefactors= colSums(data)
sizefactors
columnsums = apply(data, 2, sum)
columnsums


normalizedMatrix <- apply(data, 2, function(x)(x/sum(x))*1000000)
normalizedMatrix 
normalizedMatrixlog= log2(normalizedMatrix+1)
normalizedMatrixlog


# Assuming you have a numeric matrix named "data"
standardizedData <- matrix(NA, nrow = nrow(data), ncol = ncol(data))  # Create an empty matrix for standardized values
standardizedData
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    mean_value <- mean(data[, j])
    mean_value
    sd_value <- sd(data[, j])
    sd_value
    
    standardizedData[i, j] <- (data[i, j] - mean_value) / sd_value
    standardizedData[i, j]
  }
}


mean_value <- mean(normalizedMatrixlog)
mean_value
sd_value <- sd(normalizedMatrixlog)
sd_value
z_scores <- (normalizedMatrixlog - mean_value) / sd_value
z_scores


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library("BiocManager")
library("ComplexHeatmap")
cpm = apply(data,2,function(x)(x/sum(x))*1000000)
View(cpm)
cpm
log_transform = function(cpm){
  cpm = log2(cpm + 1)
  return(cpm)
}
log2_cpm = log_transform(cpm)
log2_cpm 

dummy = log2_cpm
dummy

gMean = apply(dummy, 1, mean)
gMean
gsd = apply(dummy, 1, sd)
gsd
calculate_zscore = function(dummy){
  z_scores <- dummy
  
  for(i in 1:nrow(dummy)){
    z_scores[i, 1:ncol(z_scores)] = (dummy[i, 1:ncol(dummy)] - gMean[i])/gsd[i]
  }
  return(z_scores)
}
m = calculate_zscore(dummy)
zsmat <- m[!rowSums(is.na(m)),]
dim(zsmat)
m
install.packages("circlize")
library("circlize")
Heatmap(zsmat[1:10, ],
col=colorRamp2(c(-2,0,2),c('orange','white','violet')))
#heatmap annotation
#choose dataset of normal vs tumore and download metafile and count data

#1.calculate variance for all genes
gene_variances <- apply(data[1:100], 1, var)
print(gene_variances)


#2.take top 100genes with most variance
top_100_genes <- names(head(sort(gene_variances, decreasing = TRUE), 100))
top_100_genes
#3.zscore those genes

top_100_expression <- zsmat[top_100_genes, ]
top_100_expression
zscored_genes <- scale(data[top_100_genes, ])
zscored_genes


zscore_expression <- scale(top_100_expression)
zscore_expression
#4. plot heatmap
install.packages('gplots')
library(gplots)   
heatmap(zscored_genes, Colv = NA, scale = "row", col = colorRampPalette(c("yellow", "white", "red"))(100),
        main = "Heatmap of Top 100 Genes (Z-Scored)", xlab = "Samples", ylab = "Genes")
