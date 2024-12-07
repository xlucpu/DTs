# ELMER Analysis: Integrative Analysis of DNA Methylation and Gene Expression
# Author: Xiaofan Lu
# Date: 2024-12-07

library(ELMER)
library(GenomicRanges)
data("Probes.motif.hg19.EPIC")

# 1. Prepare enhancer regions
# Get distal probes that are enhancers
enhancer_probes <- rownames(feature[which(feature$isEnhancer == "TRUE"),])
enhancer_probes <- intersect(enhancer_probes, rownames(annoEpic))

# Create GRanges object for enhancer regions
enhancer_granges <- annoEpic[enhancer_probes, c("CHR","MAPINFO","Strand")]
enhancer_granges$CHR <- as.character(enhancer_granges$CHR)
enhancer_granges$start <- enhancer_granges$MAPINFO
enhancer_granges$end <- enhancer_granges$start + 1
enhancer_granges$Strand <- ifelse(enhancer_granges$Strand == "F", "+", "-")

enhancer_granges <- makeGRangesFromDataFrame(
  enhancer_granges,
  keep.extra.columns = FALSE,
  ignore.strand = FALSE,
  seqinfo = NULL,
  start.field = "start",
  end.field = "end",
  strand.field = "Strand",
  starts.in.df.are.0based = FALSE
)

# 2. Process expression data
# Filter genes and prepare expression matrix
expr_mat <- tpm[intersect(rownames(tpm), Mids), rownames(annCol.meth)]
expr_mat <- expr_mat[apply(expr_mat, 1, function(x) {sum(x > 0) > 0.9*ncol(expr_mat)}),]

# Update gene annotations
Ginfo$simple <- rownames(Ginfo)
Ginfo_unique <- Ginfo[!duplicated(Ginfo$`Gene name`),]
rownames(Ginfo_unique) <- Ginfo_unique$`Gene name`

# Match genes
common_genes <- intersect(rownames(expr_mat), rownames(Ginfo_unique))
Ginfo_filtered <- Ginfo_unique[common_genes,]
expr_mat <- expr_mat[common_genes,]
rownames(expr_mat) <- Ginfo_filtered$simple # Change to ensemble ID

# 3. Create Multi-Assay Experiment object
colData <- cbind.data.frame(annCol.meth, primary = rownames(annCol.meth))
mae_enhancer <- createMAE(
  exp = expr_mat,
  met = orgmeth.desmoid[,rownames(annCol.meth)],
  save = TRUE,
  colData = colData,
  linearize.exp = FALSE,
  save.filename = file.path(res.path, "mae.enhancer.rda"),
  filter.probes = enhancer_granges,
  met.platform = "EPIC",
  genome = "hg19",
  TCGA = FALSE
)

# 4. Identify differentially methylated enhancers
sig_diff_c1 <- get.diff.meth.loose(
  data = mae_enhancer,
  group.col = "group",
  group1 = "c1",
  group2 = "c2",
  mode = "supervised",
  minSubgroupFrac = 1,
  sig.dif = 0.1,
  diff.dir = "hypo",
  cores = 1,
  dir.out = res.path,
  pvalue = 0.05,
  adjust.p = 1
)
write.table(sig_diff_c1, 
            file.path(res.path, "sig.diff.c1.enhancer.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 5. Get nearby genes for significant probes
nearby_genes <- GetNearGenes(
  data = mae_enhancer,
  probes = sig_diff_c1$probe,
  numFlankingGenes = 20
)
write.table(nearby_genes,
            file.path(res.path, "nearGenes.c1.enhancer.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 6. Identify significant probe-gene pairs
hypo_pairs <- get.pair(
  data = mae_enhancer,
  group.col = "group",
  group1 = "c1",
  group2 = "c2",
  nearGenes = nearby_genes,
  mode = "supervised",
  diff.dir = "hypo",
  permu.dir = file.path(res.path, "permu"),
  permu.size = 100000,
  raw.pvalue = 0.05,
  Pe = 0.25,
  filter.probes = TRUE,
  filter.percentage = 0.05,
  filter.portion = 0.3,
  dir.out = res.path,
  cores = 1,
  label = "hypo.c1.enhancer"
)
write.table(hypo_pairs,
            file.path(res.path, "hypo.pair.c1.enhancer.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 7. Motif enrichment analysis
enriched_motifs <- get.enriched.motif(
  data = mae_enhancer,
  probes = hypo_pairs$Probe,
  dir.out = res.path,
  probes.motif = Probes.motif.hg19.EPIC,
  background.probes = names(enhancer_granges),
  label = "hypo.c1.enhancer",
  min.incidence = 10,
  lower.OR = 1.1
)

# 8. Visualize results
pdf(file.path(fig.path, "heatmap_significant_probe_gene_pairs.pdf"),
    width = 8, height = 8)
heatmapPairs(
  data = mae_enhancer,
  group.col = "group",
  group1 = "c1",
  group2 = "c2",
  plot.distNearestTSS = FALSE,
  pairs = hypo_pairs,
  filename = NULL
)
invisible(dev.off())