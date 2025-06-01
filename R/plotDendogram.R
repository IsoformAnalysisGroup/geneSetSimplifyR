#' Plot dendogram of cluster averages
#'
#' @description A plot showing the hierarchical clustering of the cluster averages
#' (dendogram). A tool made for showing how clusters (not individual gene-sets) relate to each other.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#'
#' @return A list of ggplot objects.
#' @export
#'
#' @examples plotProportions(geneSetsList)
plotDendogram <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE
)  {
  # Check input
  if (!class(geneSetsList) == "gsList") {
    stop("'geneSetsList' should be a gsList object")
  }

  if (length(geneSetsList@cluster) == 0) {
    stop(
      "No clustering performed for 'geneSetsList'.
         Please run 'clusterGeneSets()' and 'labelClusters()' before plotting."
    )
  }

  if (length(geneSetsList@labels) == 0) {
    stop(
      "No labels created for 'geneSetsList'.
         Please run 'labelClusters()' before plotting."
    )
  }

  ### Update if not default
  if( ! is.null(resolution)) {
    geneSetsList <- changeUsedResolution(
      geneSetsList = geneSetsList,
      resolution = resolution,
      removeClusterId = removeClusterId
    )
  }

  ### Extract clusters
  df_cluster <-
    geneSetsList@cluster@active.ident %>%
    tibble::enframe() %>%
    dplyr::rename("cluster" = "value", "pathway" = "name")

  if(removeClusterId) {
    if(removeClusterId) {
      df_cluster$cluster <- stringr::str_remove(df_cluster$cluster, "^[0-9]+: ")
    }
  }
  # extract groups
  clusterVec <- unique(as.character(df_cluster$cluster))

  ### Extract distance matrix
  distMat <- 1 - geneSetsList@cluster@graphs$GSEA_snn
  distMat <- distMat[df_cluster$pathway,df_cluster$pathway]

  ### Extract average distance
  avgDistList <- lapply(
    clusterVec,
    function(group1) { # group1 <- clusterVec[1]
      lapply(
        clusterVec,
        function(group2) {  # group2 <- clusterVec[2]
          if (group1 == group2) {
            return(0)
          } else {
            avg_distance <- mean( as.vector( distMat[df_cluster$cluster == group1, df_cluster$cluster == group2] ))
            return(avg_distance)
          }
        }
    )
  })
  avgDist <- do.call(rbind, avgDistList)
  colnames(avgDist) <- clusterVec
  rownames(avgDist) <- clusterVec

  # Perform hierarchical clustering using the complete linkage method
  hcAvg <- hclust(as.dist(avgDist), method = "complete")

  # prepare for plotting
  dend <- as.dendrogram(hcAvg)
  #graphics::plot(dend, horiz=flip)

  graphics::plot(hcAvg)
}
