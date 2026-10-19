#' Run the geneSetSimplifyR workflow
#'
#' @description Wrapper function that runs the three core functions of
#' the 'geneSetSimplifyR' workflow with default parameters.
#'
#' To read more about the individual functions see their corresponding help pages:
#' \code{\link{initializeList}}, \code{\link{clusterGeneSets}}, \code{\link{labelClusters}}
#'
#' @param geneSetsList A list containing the gene sets to be analyzed.
#' The name of each element in the list should be a gene set and each element should contain
#' the genes belonging to that gene set. Alternatively, the path to a GMT file
#' with the gene sets.
#' @param geneSetsDF Optional data.frame with gene sets and enrichment scores that can be useful
#' for downstream visualization (both for effect size and direction and if multiple conditions are analyzed).
#' Should have a column named 'pathway' containing the gene sets and a column
#' named 'enrichment_score' (or 'NES') with the enrichment scores from GSEA analysis for each gene set.
#' Optionally have a column named 'pval' (or 'pvalue') with the p-value of each gene set, and
#' a column named 'source' indicating the source of each gene set. This allows for
#' analysis of multi-conditional data throughout the entire workflow and is also supported by many
#' plotting functions. See \code{\link{initializeList}} and the vignette for more info.
#' @param resolution Vector with the desired clustering resolutions to use. Higher number usually results in more clusters.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#' @param cluster_algo Clustering method applied to the SNN graph, passed to
#'   \code{clusterGeneSets()}. Either the name of a registered method (see
#'   \code{\link{listClusterMethods}}; 'louvain' by default), a
#'   \code{gsClusterMethod} from \code{\link{clusterMethod}}, or a function.
#' @param clusterArgs Named list of extra arguments passed on to the clustering
#'   method.
#' @param elbow_method Method passed to \code{clusterGeneSets()} for selecting the
#'   number of LSI dimensions. Options: 'distance'(default), 'zscore'.
#' @param K Number of nearest neighbors for the KNN graph, passed to
#'   \code{clusterGeneSets()}. `NA` (default) uses `round(sqrt(N))` for `N`
#'   gene sets.
#' @param removeFirstWord Logical. If TRUE, strip the first word (everything before
#'   the first underscore) from each gene-set name before generating cluster labels.
#'   Useful for MSigDB-style names where the first token is an author or database
#'   prefix (e.g. "GOBP", "GOMF", "HALLMARK"). Default: FALSE.
#' @param stopWords Words ignored when labeling clusters, passed to
#'   \code{labelClusters()}. Defaults to the shipped \code{\link{stopWords}}.
#' @param nLabels Number of top labels to keep per cluster, passed to
#'   \code{labelClusters()}. Default 10.
#' @return 'gsList' object with labeled clusters for all clustered resolutions.
#' @export
#'
#' @examples
#' gsListObject <- geneSetSimplifyR(
#'   geneSetsList = exampleGeneSets,
#'   geneSetsDF = exampleEnrichment,
#'   verbose = FALSE
#' )
#' gsListObject
geneSetSimplifyR <- function(
    geneSetsList,
    geneSetsDF = NULL,
    resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 4.1, 5.1),
    verbose = TRUE,
    cluster_algo = "louvain",
    clusterArgs = list(),
    elbow_method = "distance",
    K = NA,
    removeFirstWord = FALSE,
    stopWords = geneSetSimplifyR::stopWords,
    nLabels = 10
) {
  geneSetsList <- initializeList(
    geneSetsList = geneSetsList,
    geneSetsDF = geneSetsDF,
    verbose = verbose
  )

  geneSetsList <- clusterGeneSets(
    geneSetsList = geneSetsList,
    resolution = resolution,
    verbose = verbose,
    algorithm = cluster_algo,
    clusterArgs = clusterArgs,
    K = K,
    elbow_method = elbow_method
  )

  geneSetsList <- labelClusters(
    geneSetsList = geneSetsList,
    stopWords = stopWords,
    nLabels = nLabels,
    removeFirstWord = removeFirstWord,
    verbose = verbose
  )

  return(geneSetsList)
}
