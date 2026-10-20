#' @import ggplot2
#' @importFrom methods .hasSlot new
#' @importFrom stats as.dendrogram as.dist hclust mad median reorder sd
#' @importFrom utils head tail
#' @importFrom grDevices colorRampPalette
NULL

# Column names used in non-standard evaluation (dplyr / ggplot2 / tidyr).
# Declaring them here silences the "no visible binding for global variable"
# NOTE from R CMD check; behavior is unchanged.
utils::globalVariables(c(
  ".", "centrality_score", "centrality_score_scaled",
  "centrality_score_transformed", "cluster", "clusterLabel", "count",
  "enrichmentSummary", "enrichment_score", "frac", "frac_cluster",
  "gene_ids", "geneset", "genesets", "index", "label", "mean_size",
  "measure", "median_value", "n", "n_cond",
  "pathway", "rank", "scoreToSummarize", "second_label", "sortKey",
  "tableOrder", "tf_idf",
  "tf_idf_scaled", "tokenSource", "value", "weighted_score", "words",
  "x", "xend", "y", "yend"
))
