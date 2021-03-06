# Loading of additional packages
library(ggplot2)

# Definition of the working directory
setwd("")

# Definition of the data file
data.matrix <- read.table(".txt",
                          header = TRUE)

# PCA on input data
pca <- prcomp(t(data.matrix), scale = TRUE)

# Definition of data for the scree plot
pca.var <- pca$sdev ^ 2
pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
scree.data <- data.frame(Component = c(colnames(pca$x)), Percent = pca.var.per)

# Scree plot
ggplot(data = scree.data, aes(x = Component, y = Percent)) +
  geom_bar(stat = "identity", color = "Black", fill = "Black") +
  xlab("Principal Component") +
  ylab("Percentage") +
  geom_text(aes(label = paste(Percent, "%", sep = ""), vjust = -0.3)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Scree plot")

# Definition of data for the PCA plot
pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1], # PC1
                       Y = pca$x[,2]) # PC2

# Explicit names of samples for the PCA plot
pca.data["CellType"] <- c("Type 1", "Type 1", "Type 2", "Type 2")
pca.data$CellType <-as.factor(pca.data$CellType)

# PCA plot
ggplot(data = pca.data, aes(x = X, y = Y, fill = CellType)) +
  geom_point(shape = 21, size = 10, color = "Black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="right", legend.title = element_blank()) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  ggtitle("PCA plot")

# Ranking of genes contributing most to PC1
loading_scores <- pca$rotation[,1] # rotation is an output of prcomp()
gene_scores <- abs(loading_scores) # see: magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
genes_ranking <- names(gene_score_ranked["START":"END"])

# Save genes and scores
pca$rotation[genes_ranking, 1]
write.csv(pca$rotation[genes_ranking, 1], file = "genes_ranking.csv")
