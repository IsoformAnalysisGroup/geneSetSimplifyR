#' @importFrom methods .hasSlot new
#' @importFrom stats as.dendrogram as.dist hclust mad median reorder sd
#' @importFrom utils head tail
#' @importFrom grDevices colorRampPalette
NULL

# Column names used in non-standard evaluation (dplyr / ggplot2 / tidyr).
# Declaring them here silences the "no visible binding for global variable"
# NOTE from R CMD check; behaviour is unchanged.
utils::globalVariables(c(
  ".", "NES", "centrality_score", "centrality_score_scaled",
  "centrality_score_transformed", "cluster", "count", "enrichment mean",
  "enrichment_score", "frac", "frac_cluster", "gene_ids", "gene_set",
  "geneset", "genesets", "index", "label", "mean_size",
  "median_enrichment_score", "measure", "n", "org_label", "padj",
  "pathway", "tf_idf", "tf_idf_scaled", "tokenSource", "value",
  "weighted_score", "words"
))
