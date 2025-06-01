#' Run the geneSetSimplifyR workflow
#'
#' @description Wrapper function that runs the three core functions of
#' the 'geneSetSimplifyR' workflow with default parameters.
#'
#' To read more about the individual functions see their corresponding help pages:
#' \code{\link{initializeList}}, \code{\link{clusterGenesets}}, \code{\link{labelClusters}}
#'
#' @param geneSetsList A list containing the gene sets to be analyzed.
#' The name of each element in the list should be a gene set and each element should contain
#' the genes belonging to that gene set.
#' @param geneSetsDF Optional data.frame with gene sets and enrichment scores that can be usefull
#' for downstream visualization (both for effect size and direction and if multiple conditions are analyzed).
#' Should have a column named 'pathway' containing the gene sets and a column
#' named 'enrichment_score' with the enrichment scores from GSEA analysis for each gene set.
#' Optinally have a column named 'source' indicating the source of each gene set. This allows for
#' analysis of multi-condtional data throughput the entire workflow and is also supported by many
#' plotting functions. See the vignette for more info.
#' @param resolution Vector with the desired clustering resolutions to use. Higher number usually results in more clusters.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#' @return 'gsList' object with labelled clusters for all clustered resolutions.
#' @export
#'
#' @examples gsListObject <- geneSetSimplifyR(geneSetsDF, geneSetsList)
geneSetSimplifyR <- function(
    geneSetsList,
    geneSetsDF = NULL,
    resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 3.1, 4.1, 5.1),
    verbose = TRUE
) {
  geneSetsList <- initializeList(
    geneSetsList = geneSetsList,
    geneSetsDF = geneSetsDF,
    verbose = verbose
  )

  geneSetsList <- clusterGeneSets(
    geneSetsList = geneSetsList,
    resolution = resolution,
    verbose = verbose
  )

  geneSetsList <- labelClusters(
    geneSetsList = geneSetsList,
    verbose = verbose
  )

  return(geneSetsList)
}
