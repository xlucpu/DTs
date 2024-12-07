# Desmoid Tumors Analysis

## Scripts
Main R scripts for comprehensive analysis of desmoid tumors:

- `consensus_clustering.R`: Performs unsupervised consensus clustering to identify molecular subtypes in desmoid tumors
- `driver_mutation.R`: Identifies significantly mutated genes (driver mutations) using dNdScv algorithm
- `tme.R`: Analyzes tumor microenvironment (TME) by calculating TLS signatures, immune checkpoint gene expression, and TME cell population abundance
- `interprofile_correlation.R`: Calculates correlations between desmoid tumor samples and human skeletal muscle cell populations
- `ntp.R`: Performs Nearest Template Prediction (NTP) to validate subtypes in external cohorts
- `integrate_meth_expr.R`: Conducts integrative analysis of DNA methylation and gene expression using ELMER to identify methylation-regulated enhancers
- `get.diff.meth.loose.R`: A modified script of get.diff.meth() from ELMER package to allow loose stastistical criteral of using p-value less than 0.05

## Data
### Repository Data
- `Annotation.tsv`: Gene expression profiles from human extremity skeletal muscles
- `Ginfo.txt`: Gene annotation file for converting between gene symbols and Ensembl IDs

### Publication Data
The following data matrices are available in the paper's supplementary materials and Mendeley Data:
- Gene expression matrix (tpm): Available in Mendeley Data
- Mutation data (maf): Available in Table S3
- DNA methylation matrix (orgmeth.desmoid): Available in Mendeley Data
- Sample annotation (annCol): Available in Table S1
- NTP template (signature_genes): Available in Table S7

## Requirements
- R (>= 4.0.0)
- Required R packages:
  - ConsensusClusterPlus
  - dndscv
  - ELMER
  - MCPcounter
  - CMScaller
  - GSVA
  - dplyr
