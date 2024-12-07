# Calculate correlations between DT samples and human skeletal muscle cell populations
# Author: Xiaofan Lu
# Date: 2024-12-07

library(dplyr)

# Load skeletal muscle annotation data
data <- read.delim(
  file = "Annotation.tsv",
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  header = TRUE
)

# Define cell types to analyze
cell_types <- c("Adipocyte", "Fibroblast-like", "Pericyte", "SMC", "Type I", "Type II")

# Extract mean expression for each cell type
muscle_matrix <- data[, paste0(cell_types, ":mean")]

# Find common genes between muscle and tumor data
common_genes <- intersect(rownames(muscle_matrix), rownames(tpm))

# Subset both matrices to common genes
muscle_matrix <- muscle_matrix[common_genes, ]
tumor_matrix <- tpm[common_genes, rownames(annCol.merged.desmoid.expr)]

# Center expression data
muscle_matrix <- t(scale(t(muscle_matrix), scale = FALSE, center = TRUE))
tumor_matrix <- t(scale(t(tumor_matrix), scale = FALSE, center = TRUE))

# Initialize correlation matrix
cor_matrix <- matrix(
  0, 
  nrow = ncol(muscle_matrix),
  ncol = ncol(tumor_matrix), 
  dimnames = list(colnames(muscle_matrix), colnames(tumor_matrix))
) %>% as.data.frame()

# Calculate Pearson correlations between muscle types and tumor samples
for (muscle_type in rownames(cor_matrix)) {
  for (tumor_sample in colnames(cor_matrix)) {
    muscle_expr <- as.numeric(muscle_matrix[, muscle_type])
    tumor_expr <- as.numeric(tumor_matrix[, tumor_sample])
    cor_matrix[muscle_type, tumor_sample] <- cor.test(
      muscle_expr, 
      tumor_expr, 
      method = "pearson"
    )$estimate
  }
}

# cor_matrix now contains correlations between each muscle cell type and tumor sample