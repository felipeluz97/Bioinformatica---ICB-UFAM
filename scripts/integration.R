# ===============================================
# Script R: Normalização, correção de batch e PCA
# ===============================================

# Instalar
required_packages <- c("DESeq2", "edgeR", "limma", "sva", "ggplot2", "pheatmap", "RColorBrewer")
for(p in required_packages){
  if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(DESeq2)
library(edgeR)
library(limma)
library(sva)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# --------------------------
# 1. Ler dados
# --------------------------
# Suponha que temos arquivos CSV de contagem: genes x amostras
counts <- read.csv("transcriptome_counts.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1) 
# metadata deve ter colunas como: SampleID, Batch, Condition

# --------------------------
# 2. Normalização
# --------------------------
# Usando DESeq2 para normalização
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~1)  # design sem efeitos inicialmente
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized=TRUE)

# Transformação log para PCA
log_counts <- log2(norm_counts + 1)

# --------------------------
# 3. Correção de batch
# --------------------------
batch <- metadata$Batch

# 3a. Com limma::removeBatchEffect
log_counts_limma <- removeBatchEffect(log_counts, batch=batch)

# 3b. Com ComBat (sva)
log_counts_combat <- ComBat(dat=log_counts, batch=batch, par.prior=TRUE, prior.plots=FALSE)

# 3c. Usando limma com model.matrix
mod <- model.matrix(~Condition, data=metadata)
log_counts_limma_model <- removeBatchEffect(log_counts, batch=batch, design=mod)

# --------------------------
# 4. PCA
# --------------------------
plot_pca <- function(data_matrix, metadata, title){
  pca <- prcomp(t(data_matrix), scale.=TRUE)
  pca_df <- data.frame(pca$x, Batch=metadata$Batch, Condition=metadata$Condition)
  
  ggplot(pca_df, aes(x=PC1, y=PC2, color=Batch, shape=Condition)) +
    geom_point(size=3) +
    theme_minimal() +
    labs(title=title, x="PC1", y="PC2") +
    scale_color_brewer(palette="Set1")
}

# PCA antes da correção
plot_pca(log_counts, metadata, "PCA - Antes da correção de batch")

# PCA com removeBatchEffect simples
plot_pca(log_counts_limma, metadata, "PCA - Limma removeBatchEffect")

# PCA com ComBat
plot_pca(log_counts_combat, metadata, "PCA - ComBat (sva)")

# PCA com Limma + design
plot_pca(log_counts_limma_model, metadata, "PCA - Limma com design")
