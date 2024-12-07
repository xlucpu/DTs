# Script for analyzing driver mutations in desmoid tumors using dndscv
# Author: Xiaofan Lu
# Date: 2024-12-07

library(dndscv)
library(survtype)

#' Load and prepare input data
#' @param target_gene_file Path to NovoPM2.0 target gene panel file
#' @param refcds_file Path to dndscv reference CDS file
#' @return List containing target genes and reference data
prepare_input_data <- function(target_gene_file, refcds_file) {
  # Load target gene panel
  target_genes <- read.delim(
    file = target_gene_file,
    sep = "\t", 
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Load reference CDS data
  load(refcds_file)
  
  return(list(
    target_genes = setdiff(target_genes$V2, ""),
    refcds = refcds
  ))
}

#' Format mutation data for dndscv
#' @param maf Mutation data frame from Supplementary Table S3. Mutation calling by targeted gene panel sequencing for 55 desmoid tumors
#' @return Formatted mutation data frame
format_maf_for_dndscv <- function(maf) {
  data.frame(
    sampleID = maf$`ID RNA-seq`,
    chr = maf$Chr,
    pos = maf$Start,
    ref = maf$Ref,
    mut = maf$Alt,
    stringsAsFactors = FALSE
  )
}

#' Run dndscv analysis
#' @param mut_data Formatted mutation data
#' @param gene_list Target gene list
#' @return dndscv results
run_dndscv_analysis <- function(mut_data, gene_list) {
  tryCatch({
    dndscv(
      mut_data,
      cv = NULL,
      refdb = "hg19",
      gene_list = gene_list,
      outmats = TRUE,
      max_muts_per_gene_per_sample = Inf,
      max_coding_muts_per_sample = Inf
    )
  }, error = function(e) {
    stop("Error in dndscv analysis: ", e$message)
  })
}

#' Calculate mutation frequencies
#' @param maf_data Mutation data
#' @param vc_nonSyn Non-synonymous variant classifications
#' @return Data frame with mutation frequencies
calculate_mutation_frequencies <- function(maf_data, vc_nonSyn) {
  # Create mutation matrix
  mut_mat <- survtype::maf2matrix(
    as.data.frame(maf_data[which(maf_data$Variant_Classification %in% vc_nonSyn),]),
    surv.data = NULL,
    sample.name = "ID RNAseq",
    gene.name = "Hugo_Symbol",
    variant.type = "Variant_Classification"
  )
  
  # Convert to binary matrix
  mut_binary <- as.data.frame(mut_mat)
  mut_binary[mut_binary != ""] <- 1
  mut_binary[mut_binary == ""] <- 0
  mut_binary <- sapply(mut_binary, as.numeric)
  rownames(mut_binary) <- rownames(mut_mat)
  mut_binary <- as.data.frame(mut_binary)
  mut_binary$JNGR191 <- 0
  
  # Calculate frequencies
  list(
    freq = rowSums(mut_binary),
    pct = rowSums(mut_binary)/ncol(mut_binary),
    binary = mut_binary
  )
}

#' Main analysis workflow
#' @param maf Input mutation data
#' @param output_dir Output directory path
main <- function(maf, output_dir) {
  # Prepare input data
  input_data <- prepare_input_data(
    "targeted gene list 484.txt",
    "refcds_hg19.rda"
  )
  
  # Format and run dndscv
  mut_data <- format_maf_for_dndscv(maf)
  dnds_results <- run_dndscv_analysis(mut_data, input_data$target_genes)
  
  # Calculate frequencies
  mut_freqs <- calculate_mutation_frequencies(maf, vc_nonSyn)
  
  # Process results
  sel_cv <- dnds_results$sel_cv
  sel_cv <- sel_cv[sel_cv$gene_name %in% names(mut_freqs$freq),]
  rownames(sel_cv) <- sel_cv$gene_name
  
  # Add frequency information
  sel_cv$Count <- sapply(sel_cv$gene_name, function(gene) {
    nrow(maf.desmoid.filtered.rmduplicate[
      which(maf.desmoid.filtered.rmduplicate$Hugo_Symbol == gene),
    ])
  })
  sel_cv$Freq <- as.numeric(mut_freqs$freq[rownames(sel_cv)])
  sel_cv$Pct <- as.numeric(mut_freqs$pct[rownames(sel_cv)])
  sel_cv <- sel_cv[order(sel_cv$pglobal_cv),]
  
  # Write results
  write.table(
    sel_cv,
    file.path(output_dir, "desmoid_maf_hg19_dndscv_output.sig_genes_novopm_with_frequency.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

# Execute analysis
# main(maf, "results/")