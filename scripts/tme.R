# Analysis of Tumor Microenvironment in Desmoid Tumors
# Author: Xiaofan Lu
# Date: 2024-12-07

library(gsva)
library(MCPcounter)

# Define TLS signature genes
tls_signature <- c(
  # Immunoglobulin genes
  "IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", 
  "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN",
  # B cell markers 
  "CD79A", "FCRL5", "MZB1", "SSR4", "XBP1",
  # T cell markers
  "TRBC2", "IL7R",
  # Chemokines and ECM
  "CXCL12", "LUM",
  # Complement components
  "C1QA", "C7",
  # Other markers
  "CD52", "APOE", "PTLP", "PTGDS", "PIM2", "DERL3"
)

# Extract immune checkpoint gene expression
immune_checkpoint_genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3")
checkpoint_expression <- tpm[immune_checkpoint_genes, rownames(annCol)]

# Calculate TLS signature scores using ssGSEA
tls_scores <- gsva(
  ssgseaParam(
    as.matrix(tpm[, rownames(annCol)]),
    list(TLS = tls_signature)
  )
)

# Estimate abundance of TME cell populations using MCPcounter
tme_populations <- MCPcounter.estimate(
  tpm,
  featuresType = "HUGO_symbols"
)

# Results now stored in:
# - checkpoint_expression: Expression of immune checkpoint genes
# - tls_scores: TLS signature scores
# - tme_populations: Abundance of TME cell populations