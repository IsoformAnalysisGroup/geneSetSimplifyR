#' Plot label information for a particular cluster
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#' @param resolution Float. The clustering resolution to use as the default.
#' The labels of this resolution will be set as the 'active.idents' in the Seurat object containing the clustering,
#' and thus this is the resolution used for the different plotting functions.
#' The 'resolution' must have been used as a resolution in the 'clusterGeneSets()' function.
#' @param clusterNumber Which cluster to plot. Here the cluster number is used
#' @param newClusterName The new name of the cluster
#'
#' @return 'gsList' object with a single cluster label updated. The original names are also stored.
#' @export
#'
#' @examples
#' exampleGsList <- updateClusterLabels(geneSetsList = exampleGsList, 
#' clusterNumber = 0, newClusterName = 'cancer response')
updateClusterLabels <- function(
    geneSetsList,
    resolution = NULL,
    clusterNumber,
    newClusterName
) {
  if( ! is.null(resolution)) {
    geneSetsList <- changeUsedResolution(
      geneSetsList = geneSetsList,
      resolution = resolution,
      removeClusterId = FALSE
    )
  }


  if( ! is.null(resolution)) {
    resolutionColumn <-
      paste0("GSEA_snn_res.", resolution)
  } else {
    resolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]
  }
  ### Save original
  geneSetsList@labels[[resolutionColumn]]$labelsDefault$org_label <-
    geneSetsList@labels[[resolutionColumn]]$labelsDefault$label

  ### Overwrite with new
  toModify <- which(
    geneSetsList@labels[[resolutionColumn]]$labelsDefault$cluster == clusterNumber
  )[1]

  geneSetsList@labels[[resolutionColumn]]$labelsDefault$label[toModify] <- newClusterName

  ### Update stored labes
  defaultLabels <-
    geneSetsList@labels[[resolutionColumn]][["labelsDefault"]] %>%
    dplyr::mutate(cluster = as.numeric(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(weighted_score, n = 1, with_ties = F) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = paste0(cluster, ": ", label)) %>%
    dplyr::pull(label)

  currentClusterLabels <-
    levels(geneSetsList@cluster@active.ident)

  geneSetsList@cluster@active.ident <-
    plyr::mapvalues(x = geneSetsList@cluster@active.ident,
                    from = currentClusterLabels,
                    to = defaultLabels)

  return(geneSetsList)
}

