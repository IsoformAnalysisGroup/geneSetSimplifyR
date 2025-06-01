#' Plot cluster tree
#'
#' @description Graphs a clustering tree plot showing t
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param node_alpha Numeric value indicating the alpha of nodes
#' @param node_text_sixe Numeric value indicating node text size
#' @param node_size_range Numeric vector with two values indicating the
#' minimum and maximum point size of nodes, respectively
#' @param edge_width Numeric value indicating the width of plotted edges
#' @param node_text_angle Numeric value indicating the rotation of node labels
#' @param resolutionsToPlot Resolutions to include in plot. Use Null to plot all. Default is NULL
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples plotClustertree(geneSetsList)
plotClustree <- function(
    geneSetsList,
    node_alpha = 1,
    node_text_size = 2.5,
    node_size_range = c(4, 20),
    edge_width = 0.7,
    node_text_angle = 30,
    resolutionsToPlot = NULL
) {
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

  resolutions <-
    grep("GSEA_snn_res.",
         names(geneSetsList@cluster@meta.data),
         value = TRUE)

  ### Subset to resolutions
  if(!is.null(resolutionsToPlot)) {
    resolutionsToUse <- paste0(
      "GSEA_snn_res.",
      resolutionsToPlot
    )

    if(! all(resolutionsToUse %in% resolutions)) {
      stop('All \'resolutionsToUse\' must be among the resolutions analyzed')
    }

    notToAnalyze <- setdiff(resolutions, resolutionsToUse)

    resolutions <- intersect(resolutions, resolutionsToUse)

    ### Remove from seurat object
    geneSetsList@cluster@meta.data[,notToAnalyze] <-list(NULL)

  }

  for (resolution in resolutions) {
    labels <-
      geneSetsList@labels[[resolution]][["labelsDefault"]] %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(weighted_score, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      tidyr::unite("label", cluster:label, sep = ": ") %>%
      dplyr::select(label) %>%
      dplyr::pull()

    levels(geneSetsList@cluster@meta.data[[resolution]]) <-
      labels
  }
  p <- clustree::clustree(
    geneSetsList@cluster,
    prefix = "GSEA_snn_res.",
    label_nodes = TRUE,
    node_alpha = node_alpha,
    node_text_size = node_text_size,
    node_size_range = node_size_range,
    edge_width = edge_width,
    node_text_angle = node_text_angle
  )

  return(p)
}
