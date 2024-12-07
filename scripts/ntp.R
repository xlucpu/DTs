# Nearest Template Prediction (NTP) Analysis for DT Subtype Classification
# Author: Xiaofan Lu
# Date: 2024-12-07

library(CMScaller)

# Load subtype signature genes (top 200 upregulated genes per subtype)
# not run
# signature_genes <- read.table(
#   "Table S7. Top 200 significantly upregulated genes in each subtype"
# )

# Format template for NTP
templates <- signature_genes[, c("Gene", "Subtype")]
colnames(templates) <- c("probe", "class")

# Scale expression data
scaled_expression <- t(scale(t(expr.validation), scale = TRUE, center = TRUE))

# Perform NTP classification
subtype_predictions <- CMScaller::ntp(
  scaled_expression,
  templates = templates,
  doPlot = TRUE,
  seed = 20000112
)

# subtype_predictions now contains predicted subtype labels for validation samples