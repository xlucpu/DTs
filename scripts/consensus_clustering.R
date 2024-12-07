# Consensus clustering analysis for desmoid tumor transcriptomic subtypes
# Author: Xiaofan Lu
# Date: 2024-12-07

library(ConsensusClusterPlus)
library(pheatmap)
library(NMF)

#' Preprocess expression data
#' @param tpm TPM expression matrix
#' @param annCol Clinical annotation data from Supplementary Table S1. Clinical and molecular characteristics of 76 desmoid tumors analyzed in this study
#' @param n_features Number of top variable features to select
#' @return Preprocessed expression matrix
preprocess_expression_data <- function(tpm, annCol, n_features = 1000) {
  # Align data with annotation
  indata <- tpm[, rownames(annCol)]
  
  # Remove NA values and filter low-expressed genes
  indata <- na.omit(indata)
  indata <- indata[apply(indata, 1, function(x) {
    sum(x > 0) > 0.9 * ncol(indata)
  }), ]
  
  # Select top variable genes
  var.sd <- apply(indata, 1, mad)
  var.sel <- names(sort(var.sd, decreasing = TRUE)[1:n_features])
  
  # Center the data
  indata <- t(scale(t(indata[var.sel,]), scale = FALSE, center = TRUE))
  
  return(indata)
}

#' Run consensus clustering
#' @param indata Preprocessed expression matrix
#' @param max_k Maximum number of clusters to try
#' @param reps Number of subsampling repetitions
#' @param output_dir Output directory for plots
#' @return Consensus clustering results
run_consensus_clustering <- function(indata, max_k = 6, reps = 500, output_dir) {
  cc.mrna <- ConsensusClusterPlus(
    d = as.matrix(indata),
    maxK = max_k,
    reps = reps,
    pItem = 1,
    pFeature = 0.9,
    clusterAlg = "pam",
    innerLinkage = "ward.D",
    finalLinkage = "ward.D",
    distance = "pearson",
    seed = 20000112,
    title = file.path(output_dir, "consensusmRNA"),
    plot = "pdf"
  )
  
  return(cc.mrna)
}

#' Assign cluster labels and update annotation
#' @param cc_results Consensus clustering results
#' @param k Number of clusters to use
#' @param annCol Clinical annotation data from Supplementary Table S1. Clinical and molecular characteristics of 76 desmoid tumors analyzed in this study
#' @return Updated annotation with cluster labels
assign_clusters <- function(cc_results, k, annCol) {
  group <- cc_results[[k]]$consensusClass
  group <- paste0("C", group)
  names(group) <- colnames(indata)
  annCol$Subtype <- group[rownames(annCol)]
  return(annCol)
}

#' Plot consensus matrix heatmap
#' @param cc_results Consensus clustering results
#' @param k Number of clusters
#' @param annCol Updated annotation with cluster labels
#' @param cellwidth Width of heatmap cells
plot_consensus_matrix <- function(cc_results, k, annCol, cellwidth = 10) {
  # Prepare plotting data
  hcs <- cc_results[[k]]$consensusTree
  plotdata <- cc_results[[k]]$consensusMatrix
  dimnames(plotdata) <- list(colnames(indata), colnames(indata))
  
  # Generate heatmap
  pheatmap(
    mat = as.matrix(plotdata),
    border_color = NA,
    color = NMF:::ccRamp(heatmap.BlWtRd2, 64),
    show_rownames = FALSE,
    show_colnames = TRUE,
    cluster_rows = hcs,
    cluster_cols = hcs,
    annotation_col = annCol,
    name = "mRNA",
    cellwidth = cellwidth,
    cellheight = 100/nrow(plotdata)
  )
}

#' Main analysis workflow
#' @param tpm TPM expression matrix
#' @param annCol Clinical annotation data
#' @param output_dir Output directory
#' @param n_features Number of top variable features
#' @param k Number of clusters
main <- function(tpm, annCol, output_dir, n_features = 1000, k = 5) {
  # Preprocess data
  indata <- preprocess_expression_data(tpm, annCol, n_features)
  
  # Run consensus clustering
  cc_results <- run_consensus_clustering(indata, output_dir = output_dir)
  
  # Assign clusters and update annotation
  annCol <- assign_clusters(cc_results, k, annCol)
  
  # Plot results
  plot_consensus_matrix(cc_results, k, annCol)
  
  return(list(
    clusters = annCol$Subtype,
    cc_results = cc_results
  ))
}

# Execute analysis
# results <- main(tpm, annCol, "results/")