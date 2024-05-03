## loading packages

library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(ggthemes)
library(tidyverse)
library(pheatmap)

## reading data
data <- read.delim("RNA2_vs_noguide.txt",sep="\t",header=T)

data2 <- data[,c(2,18:23,25,26)]

data2 <- unique(data2)

data3 <- data2[which(data2$Gene.Symbol!=""),]

data4 <- data3[!duplicated(data3$Gene.Symbol),]

data4$miss1 <- apply(data4[,2:4],1,function(x) {sum(is.na(x))})
data4$miss2 <- apply(data4[,5:7],1,function(x) {sum(is.na(x))})

data4$variance <- apply(data4[,2:7],1,var,na.rm=TRUE)

## select genes with top 100 variances and showing values at at least two samples for each condition
data5 <- head(arrange(data4[which(data4$miss1<=1 & data4$miss2<=1),],desc(variance)), n = 100)


numerical_data <- data5[,2:7]
names(numerical_data) <- gsub("Abundances..Normalized...","",names(numerical_data))
rownames(numerical_data) <- data5$Gene.Symbol


# Calculate distances between samples
sampleDists <- dist(t(data4[,2:7]))

# Plot inter-sample distances
old.par <- par(no.readonly=T)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- gsub("\\.\\.",".",gsub("Abundances..Normalized...","",rownames(sampleDistMatrix)))
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p.hm <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 color = colorRampPalette(c("red", "white", "blue"))(100))

pdf("RNA2_distance_heat_map.pdf", p.hm,width=7,height=6)
p.hm
dev.off()

p.hm


## Compute the correlation matrix
corr_matrix <- cor(numerical_data, use="complete.obs")
ggcorrplot(corr_matrix,hc.order = TRUE)

## Applying PCA
data.pca <- princomp(corr_matrix)
summary(data.pca)


## Visualization of the principal components 
# Scree Plot
fviz_eig(data.pca, addlabels = TRUE)

## Biplot of the attributes
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")


## Contribution of each variable 
fviz_cos2(data.pca, choice = "var", axes = 1:2)

## Biplot combined with cos2 
p <- fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)


pdf("RNA2_pca.pdf",height=6,width=6)
p
dev.off()


## ploting PCA with ggplot2

df <- as.data.frame(data.pca$scores)
df$class <- gsub("^.*\\.\\.","",rownames(df))


PoV <- data.pca$sdev^2/sum(data.pca$sdev^2)

PC1 <- paste(round(100*PoV[1], 2), "%", sep="")
PC2 <- paste(round(100*PoV[2], 2), "%", sep="")

rownames(df) <- gsub("\\.\\.",".",rownames(df))

p2 <- ggplot(df) +
  aes(Comp.1, Comp.2, color = class, shape = class, label=rownames(df)) + 
  geom_point(size = 2) + 
  geom_text(nudge_x = 0, nudge_y=-0.005) +
  #coord_fixed() +
  xlab(paste0("PC1: ",PC1))+ 
  ylab(paste0("PC2: ",PC2))+
  xlim(limits=c(-0.12,0.16))+
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  theme_bw()

p2

pdf("RNA2_pca2.pdf",height=4.5,width=6)
p2
dev.off()

