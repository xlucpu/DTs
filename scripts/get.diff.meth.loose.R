library(plyr)
get.diff.meth.loose <- function (data, diff.dir = "hypo", cores = 1, mode = "unsupervised", 
          minSubgroupFrac = 0.2, pvalue = 0.05, adjust.p = 0.05, group.col, min.samples = 5, 
          group1, group2, test = t.test, sig.dif = 0.3, dir.out = "./", 
          save = TRUE) 
{
  if (is.null(getMet(data))) 
    stop("Cannot identify differential DNA methylation region without DNA methylation data.")
  if (nrow(colData(data)) == 0) {
    stop("Sample information data to do differential analysis.")
  }
  else if (missing(group.col)) {
    stop("Please colData.col should be specified, labeling two group of sample for comparison. See colnames(colData(data)) for possibilities")
  }
  else if (!group.col %in% colnames(colData(data))) {
    stop("Group column not found in phenotypic data and meta-data of the object. See values with colData(data)")
  }
  else if (missing(group1) | missing(group2)) {
    if (length(unique(colData(data)[, group.col])) < 2) {
      stop("Group column should have at least 2 distinct group labels for comparison.")
    }
    else if (length(unique(colData(data)[, group.col])) > 
             2) {
      stop("Please your object must have only two groups. We found more than two and this might impact the next analysis steps.")
    }
    else {
      groups <- colData(data)[, group.col]
      group1 <- unique(groups)[1]
      group2 <- unique(groups)[2]
      message(paste0("Group 1: ", group1, "\nGroup 2: ", 
                     group2))
    }
  }
  else if (!group1 %in% unique(colData(data)[, group.col])) {
    stop(group1, " not found in ", group.col)
  }
  else if (!group2 %in% unique(colData(data)[, group.col])) {
    stop(group2, " not found in ", group.col)
  }
  if (!diff.dir %in% c("hypo", "hyper", "both")) 
    stop("diff.dir optiosn are hypo, hyper or both")
  if (diff.dir %in% c("both")) 
    diff.dir <- NA
  parallel <- FALSE
  if (cores > 1) {
    if (cores > detectCores()) 
      cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  Top.m <- ifelse(diff.dir == "hyper", TRUE, FALSE)
  if (is.na(Top.m) & minSubgroupFrac < 1) {
    message("Two tailed test should be performed with all samples")
    minSubgroupFrac <- 1
  }
  if (mode == "supervised" & minSubgroupFrac < 1) {
    message("Supervised mode will use all samples from boths groups. Setting argument minSubgroupFrac to 1")
    minSubgroupFrac <- 1
  }
  counts <- plyr::count(MultiAssayExperiment::colData(data)[, 
                                                            group.col])
  message(paste0("ELMER will search for probes ", ifelse(is.na(diff.dir), 
                                                         "differently ", diff.dir), "methylated in group ", group1, 
                 " (n:", subset(counts, counts$x == group1)$freq, ")", 
                 " compared to ", group2, " (n:", subset(counts, counts$x == 
                                                           group2)$freq, ")"))
  message(paste0("ooo Arguments ooo"))
  message(paste0("o Number of probes: ", nrow(getMet(data))))
  message(paste0("o Beta value difference cut-off: ", sig.dif))
  message(paste0("o P-value cut-off: ", pvalue))
  message(paste0("o FDR cut-off: ", adjust.p))
  message(paste0("o Mode: ", mode))
  message(paste0("o % of samples per group in each comparison: ", 
                 minSubgroupFrac))
  message(paste0("o Min number of samples per group in each comparison: ", 
                 min.samples))
  message(paste0("o Nb of samples group1 in each comparison: ", 
                 ifelse(round(subset(counts, counts$x == group1)$freq * 
                                minSubgroupFrac) > min.samples, round(subset(counts, 
                                                                             counts$x == group1)$freq * minSubgroupFrac), min(min.samples, 
                                                                                                                              subset(counts, counts$x == group1)$freq))))
  message(paste0("o Nb of samples group2 in each comparison: ", 
                 ifelse(round(subset(counts, counts$x == group2)$freq * 
                                minSubgroupFrac) > min.samples, round(subset(counts, 
                                                                             counts$x == group2)$freq * minSubgroupFrac), min(min.samples, 
                                                                                                                              subset(counts, counts$x == group2)$freq))))
  message(paste0("Output direction: ", dir.out))
  message(paste0("ooooooooooooooooo"))
  groups.info <- colData(data)[getMetSamples(data), group.col]
  met <- assay(getMet(data))
  probes <- rownames(met)
  out <- alply(.data = met, .margins = 1, .fun = function(x) {
    ELMER:::Stat.diff.meth(percentage = minSubgroupFrac, meth = x, 
                   min.samples = min.samples, groups = groups.info, 
                   group1 = group1, test = test, group2 = group2, Top.m = Top.m)
  }, .progress = "time", .parallel = parallel, .paropts = list(.errorhandling = "pass"))
  out <- do.call(rbind, out)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$probe <- probes
  diffCol <- paste0(gsub("[[:punct:]]| ", ".", group1), "_Minus_", 
                    gsub("[[:punct:]]| ", ".", group2))
  out$adjust.p <- p.adjust(as.numeric(out$PP), method = "BH")
  out <- out[, c("probe", "PP", "MeanDiff", "adjust.p")]
  colnames(out) <- c("probe", "pvalue", diffCol, "adjust.p")
  rownames(out) <- out$probe
  if (save) {
    message("Saving results")
    dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
    ylab <- ifelse(is.na(diff.dir), " (FDR corrected P-values) [two tailed test]", 
                   " (FDR corrected P-values) [one tailed test]")
    TCGAVisualize_volcano(x = as.data.frame(out)[, grep("Minus", 
                                                        colnames(out), value = T)], y = out$adjust.p, title = paste0("Volcano plot - Probes ", 
                                                                                                                     ifelse(is.na(diff.dir), "differently ", diff.dir), 
                                                                                                                     "methylated in ", group1, " vs ", group2, "\n"), 
                          filename = sprintf("%s/volcanoPlot.probes.%s.png", 
                                             dir.out, ifelse(is.na(diff.dir), "two_tailed", 
                                                             diff.dir)), label = c("Not Significant", paste0("Hypermethylated in ", 
                                                                                                             group1), paste0("Hypomethylated in ", group1)), 
                          ylab = bquote(-Log[10] ~ .(ylab)), xlab = expression(paste("DNA Methylation difference (", 
                                                                                     beta, "-values)")), x.cut = sig.dif, y.cut = pvalue)
    write_csv(x = out, file = sprintf("%s/getMethdiff.%s.probes.csv", 
                                      dir.out, ifelse(is.na(diff.dir), "both", diff.dir)))
    write_csv(x = out[out$pvalue < pvalue & out$adjust.p < adjust.p & abs(out[, diffCol]) > 
                        sig.dif & !is.na(out$pvalue), ], file = sprintf("%s/getMethdiff.%s.probes.significant.csv", 
                                                                          dir.out, ifelse(is.na(diff.dir), "both", diff.dir)))
  }
  result <- out[out$pvalue < pvalue & out$adjust.p < adjust.p & abs(out[, diffCol]) > 
                  sig.dif & !is.na(out$pvalue), ]
  if (nrow(result) == 0) {
    message("No relevant probes found")
  }
  else {
    message(paste0("Number of relevant probes found: ", nrow(result)))
  }
  return(result)
}