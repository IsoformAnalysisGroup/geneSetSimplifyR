#' Rename a cluster
#'
#' @description Replaces the label of one cluster with a name of your choice.
#'   The original label is kept and reported as `initial_cluster_label` by
#'   `extractTable(summarize = FALSE)`, also after several renames.
#'
#' @param geneSetsList 'gsList' object that has been clustered and labeled.
#' @param resolution Numeric. The clustering resolution of the cluster to
#'   rename. It must have been used as a resolution in `clusterGeneSets()`.
#'   When given, this resolution also becomes the active resolution: plotting
#'   and table functions called with `resolution = NULL` afterwards use it.
#'   Default `NULL` renames a cluster of the active resolution.
#' @param clusterNumber The number of the cluster to rename, as shown in front
#'   of its label (e.g. `0` for "0: Immune Response").
#' @param newClusterName The new name of the cluster.
#'
#' @return 'gsList' object with the cluster label updated.
#' @export
#'
#' @examples
#' gsListObject <- updateClusterLabels(
#'   geneSetsList = exampleGsList,
#'   clusterNumber = 0,
#'   newClusterName = 'cancer response'
#' )
#' head(extractTable(gsListObject, summarize = FALSE))
updateClusterLabels <- function(
    geneSetsList,
    resolution = NULL,
    clusterNumber,
    newClusterName
) {
  checkGsList(geneSetsList)
  geneSetsList <- upgradeGsList(geneSetsList)

  if (missing(clusterNumber) || length(clusterNumber) != 1L ||
      is.na(clusterNumber)) {
    stop("'clusterNumber' should be a single cluster number.")
  }
  if (missing(newClusterName) || !is.character(newClusterName) ||
      length(newClusterName) != 1L || is.na(newClusterName) ||
      !nzchar(trimws(newClusterName))) {
    stop("'newClusterName' should be a single non-empty string.")
  }

  resolutionColumn <- resolveResolutionColumn(geneSetsList, resolution)

  clusterIds <- levels(geneSetsList@cluster@meta.data[[resolutionColumn]])
  if (!as.character(clusterNumber) %in% clusterIds) {
    stop(
      "Cluster ", clusterNumber, " does not exist at resolution ",
      sub("^GSEA_snn_res\\.", "", resolutionColumn), ". Available clusters: ",
      paste(clusterIds, collapse = ", "), "."
    )
  }

  labelsDefault <- geneSetsList@labels[[resolutionColumn]]$labelsDefault

  ### Save the original labels the first time a label is changed
  if (!'org_label' %in% colnames(labelsDefault)) {
    labelsDefault$org_label <- labelsDefault$label
  }

  ### Warn if the new name is already used by another cluster
  otherLabels <- alignedClusterLabels(
    geneSetsList = geneSetsList,
    resolutionColumn = resolutionColumn,
    clusterIds = setdiff(clusterIds, as.character(clusterNumber)),
    removeClusterId = TRUE
  )
  if (newClusterName %in% otherLabels) {
    warning(
      "Another cluster is already labeled '", newClusterName, "'. ",
      "The two clusters stay separate but can be hard to tell apart.",
      call. = FALSE
    )
  }

  ### Overwrite the top label of the cluster (the one that is displayed)
  inCluster <- which(
    as.character(labelsDefault$cluster) == as.character(clusterNumber)
  )
  if (length(inCluster)) {
    toModify <- inCluster[which.max(labelsDefault$weighted_score[inCluster])]
    labelsDefault$label[toModify] <- newClusterName
  } else {
    ### Clusters without any label get one
    newRow <- labelsDefault[0, , drop = FALSE]
    newRow[1, 'label'] <- newClusterName
    newRow[1, 'cluster'] <- as.numeric(clusterNumber)
    newRow[1, 'weighted_score'] <- Inf
    labelsDefault <- dplyr::bind_rows(labelsDefault, newRow)
  }

  geneSetsList@labels[[resolutionColumn]]$labelsDefault <- labelsDefault

  ### Make this resolution the active one and refresh the identities
  geneSetsList@metadata$active_resolution <- resolutionColumn
  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE
  )

  return(geneSetsList)
}
