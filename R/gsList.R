#' @keywords internal
"_PACKAGE"

#' The gsList class
#'
#' @description A `gsList` holds the input and results of the geneSetSimplifyR
#'   workflow. It is created by [initializeList()], filled by
#'   [clusterGeneSets()] and [labelClusters()], or all at once by
#'   [geneSetSimplifyR()]. Use [extractTable()] to get the results as a table.
#'
#' @slot genesetsDataframe A data.frame with genes as rows and gene sets as
#'   columns; 1 if the gene is in the gene set, else 0.
#' @slot enrichmentScore A data.frame with one row per gene set and source:
#'   `pathway`, `source`, `enrichment_score`, `pval`, `nlog10pval` and any
#'   other column supplied to [initializeList()].
#' @slot cluster A Seurat object with the LSI, UMAP, SNN graph and the cluster
#'   of every gene set for each resolution (`GSEA_snn_res.<resolution>`
#'   metadata columns).
#' @slot labels A list with, for each resolution, the cluster labels
#'   (`labelsDefault`), the TF-IDF scores of the label candidates
#'   (`labelsTFIDF`) and the most central gene sets (`labelsCentrality`).
#' @slot metadata A list with parameters and results of the workflow, e.g.
#'   the number of nearest neighbors, the chosen resolution and silhouette
#'   scores.
#'
#' @return `gsList` objects are created by [initializeList()], not by calling
#'   the class generator directly. `show()` prints a summary of the object and
#'   returns it invisibly.
#'
#' @exportClass gsList
#' @examples
#' exampleGsList
#' names(exampleGsList@metadata)
gsList <- setClass(
  "gsList",
  slots = c(
    genesetsDataframe = "data.frame",
    enrichmentScore = "data.frame",
    cluster = "ANY",
    labels = "list",
    metadata = "list"
  ),
  prototype = list(
    metadata = list()
  )
)

#' Add slots missing from objects created by older package versions
#'
#' @param geneSetsList A 'gsList' object.
#' @return The object, with an empty `metadata` slot added if it had none.
#' @keywords internal
upgradeGsList <- function(geneSetsList) {
  if (.hasSlot(geneSetsList, "metadata")) {
    return(geneSetsList)
  }
  methods::new(
    "gsList",
    genesetsDataframe = geneSetsList@genesetsDataframe,
    enrichmentScore = geneSetsList@enrichmentScore,
    cluster = geneSetsList@cluster,
    labels = geneSetsList@labels,
    metadata = list()
  )
}
