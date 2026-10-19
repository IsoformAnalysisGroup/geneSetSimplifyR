#' Plot cluster tree
#'
#' @description Graphs a clustering tree showing how the clusters of
#'   consecutive resolutions relate to each other: each node is a cluster
#'   (labeled with its top label), and edges show how many gene sets move
#'   from a cluster at one resolution to a cluster at the next.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param node_alpha Numeric value indicating the alpha of nodes
#' @param node_text_size Numeric value indicating node text size
#' @param node_size_range Numeric vector with two values indicating the
#' minimum and maximum point size of nodes, respectively
#' @param edge_width Numeric value indicating the width of plotted edges
#' @param node_text_angle Numeric value indicating the rotation of node labels
#' @param resolutionsToPlot Resolutions to include in plot. Use NULL to plot all. Default is NULL
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' plotClustree(exampleGsList)
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
  checkGsList(geneSetsList)

  resolutions <- availableResolutionColumns(geneSetsList)

  ### Subset to resolutions
  if(!is.null(resolutionsToPlot)) {
    resolutionsToUse <- resolutionColumnName(resolutionsToPlot)

    if(! all(resolutionsToUse %in% resolutions)) {
      stop(
        "All 'resolutionsToPlot' must be among the resolutions analyzed: ",
        paste(sub("^GSEA_snn_res\\.", "", resolutions), collapse = ", "), "."
      )
    }

    notToAnalyze <- setdiff(resolutions, resolutionsToUse)

    resolutions <- intersect(resolutions, resolutionsToUse)

    ### Remove from seurat object
    geneSetsList@cluster@meta.data[,notToAnalyze] <-list(NULL)

  }

  for (resolution in resolutions) {
    ### Top label per cluster, keyed by its cluster id
    labelMap <-
      geneSetsList@labels[[resolution]][["labelsDefault"]] %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(weighted_score, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(label = paste0(cluster, ": ", label)) %>%
      dplyr::select(cluster, label)

    ### Match labels to the factor's levels BY cluster id, not by position.
    ### The label vector's order (e.g. numeric group order) need not match the
    ### factor's level order, so a positional `levels<-` can attach a label to
    ### the wrong cluster's node. Matching on the id keeps them aligned.
    clusterFactor <- geneSetsList@cluster@meta.data[[resolution]]
    matchedLabels <- labelMap$label[
      match(as.character(levels(clusterFactor)),
            as.character(labelMap$cluster))
    ]

    ### Fall back to the raw cluster id for any level without a label
    missingLabel <- is.na(matchedLabels)
    matchedLabels[missingLabel] <- levels(clusterFactor)[missingLabel]

    levels(clusterFactor) <- matchedLabels
    geneSetsList@cluster@meta.data[[resolution]] <- clusterFactor
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

  # clustree's edge color scale names its guide ("edge_colourbar"), which
  # ggplot2 can only find when ggraph is attached. Use the guide object itself
  # so the plot prints without changing the user's search path.
  for (scale in p$scales$scales) {
    if (any(grepl("^edge_colour", scale$aesthetics)) &&
        is.character(scale$guide) && scale$guide == "edge_colourbar") {
      scale$guide <- ggraph::guide_edge_colourbar()
    }
  }

  return(p)
}
