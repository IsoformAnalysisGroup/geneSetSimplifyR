#' Scatter plot of UMAP dimensionality reduction
#'
#' @description Graphs the output of the UMAP dimensionality reduction as a scatterplot
#' using the 'DimPlot()' function from the Seurat package.
#' Gene sets are colored and labeled by their cluster.
#' When there is only one source in the data (default),
#' each point will be shaped by the direction of the enrichment score of the gene set.
#' When there are two sources in the data,
#' each point will be shaped by the source indicating whether the gene set was found
#' significant in either or both sources. With three or more sources, points
#' are not shaped.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param interactive Logical, indicating whether to plot an interactive UMAP plot. Default: FALSE.
#' If TRUE, an interactive plot is made with 'HoverLocator()' from the Seurat
#' package; hovering over a point shows its gene set and cluster.
#' @param highlightClusterNo A vector with the cluster numbers of the clusters (at the resolution) to highlight in the plot.
#' @param pointSize Size of points. Default is 2.
#' @param pointAlpha Transparency of points. Default is 1/3.
#' @param labelSize Size of the cluster labels. Default is 3.
#' @param labelAlpha Transparency of the cluster labels. Default is 1.
#' @param highlightPointAlpha Transparency of the points of highlighted clusters. Default is 2/3.
#' @param highlightPointSize Size of the points of highlighted clusters. Default is 4.
#' @param highlightLabelAlpha Transparency of the labels of highlighted clusters. Default is 1.
#' @param highlightLabelSize Size of the labels of highlighted clusters. Default is 5.
#'
#' @return A ggplot object (a patchwork from 'Seurat::DimPlot()'), or a
#'   plotly object when `interactive = TRUE`.
#' @export
#'
#' @examples
#' plotUMAP(exampleGsList)
#' plotUMAP(exampleGsList, highlightClusterNo = 0, removeClusterId = TRUE)
plotUMAP <- function(
    ### High level
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    interactive = FALSE,
    highlightClusterNo = NULL,
    ### Detailed
    pointSize = 2,
    pointAlpha = 1/3,
    labelSize = 3,
    labelAlpha = 1,
    highlightPointAlpha = 2/3,
    highlightPointSize = 4,
    highlightLabelAlpha = 1,
    highlightLabelSize = 5
) {
  # Check input
  checkGsList(geneSetsList)

  clusterNameToHighlight <- highlightedClusterNames(
    geneSetsList = geneSetsList,
    resolution = resolution,
    highlightClusterNo = highlightClusterNo,
    removeClusterId = removeClusterId
  )

  ### Set the cluster resolution and labels to plot
  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = resolution,
    removeClusterId = removeClusterId
  )

  ### Shape the points by source or by direction of enrichment
  nCond <- countConditions(geneSetsList)
  geneSetIds <- geneSetsList@cluster@meta.data$gene_set
  shapeBy <- NULL
  shapeTitle <- NULL
  hoverVars <- c("ident", "nCount_GSEA")

  if( nCond == '2' ) {
    sourcesPerGeneSet <- table(geneSetsList@enrichmentScore$pathway)
    firstSource <- geneSetsList@enrichmentScore$source[match(
      geneSetIds, geneSetsList@enrichmentScore$pathway
    )]
    geneSetSource <- ifelse(
      as.vector(sourcesPerGeneSet[geneSetIds]) > 1,
      'Both',
      firstSource
    )

    geneSetsList@cluster <- SeuratObject::AddMetaData(
      object = geneSetsList@cluster,
      metadata = stats::setNames(geneSetSource, colnames(geneSetsList@cluster)),
      col.name = 'source'
    )
    shapeBy <- 'source'
    shapeTitle <- 'Gene set source'
    hoverVars <- c("ident", "source", "nCount_GSEA")
  }

  hasEnrichmentScore <-
    ! all(is.na( geneSetsList@enrichmentScore$enrichment_score ))

  if( nCond == '1' && hasEnrichmentScore ) {
    # Matched by gene-set name; the rows of 'enrichmentScore' need not be in
    # the same order as the gene sets in the Seurat object.
    enrichmentScore <- geneSetsList@enrichmentScore$enrichment_score[match(
      geneSetIds, geneSetsList@enrichmentScore$pathway
    )]
    enrichmentScoreSign <- dplyr::case_when(
      enrichmentScore > 0 ~ 'Positive',
      enrichmentScore < 0 ~ 'Negative',
      enrichmentScore == 0 ~ 'Zero',
      TRUE ~ 'Missing'
    )

    geneSetsList@cluster <- Seurat::AddMetaData(
      object = geneSetsList@cluster,
      metadata = data.frame(
        enrichment_score_sign = enrichmentScoreSign,
        enrichment_score = enrichmentScore,
        row.names = colnames(geneSetsList@cluster)
      )
    )
    shapeBy <- 'enrichment_score_sign'
    shapeTitle <- 'Enrichment score'
    hoverVars <- c("ident", "enrichment_score", "nCount_GSEA")
  }

  ### Make plot
  pUMAP <- Seurat::DimPlot(
    geneSetsList@cluster,
    label = TRUE,
    label.box = TRUE,
    repel = TRUE,
    label.size = labelSize,
    pt.size = pointSize,
    shape.by = shapeBy,
    group.by = 'ident'
  ) +
    theme_bw() +
    guides(colour = "none") +
    labs(title = NULL, x = 'UMAP1', y = 'UMAP2')

  if (!is.null(shapeBy)) {
    pUMAP <- pUMAP +
      theme(legend.position = "bottom") +
      labs(shape = shapeTitle)
  }

  ### Set transparency and size
  pointIsHighlighted <- pUMAP[[1]]$data$ident %in% clusterNameToHighlight
  labelIsHighlighted <- pUMAP[[1]]$layers[[2]]$data$ident %in% clusterNameToHighlight

  # points
  pUMAP[[1]]$layers[[1]]$aes_params$alpha <-
    ifelse(pointIsHighlighted, highlightPointAlpha, pointAlpha)
  pUMAP[[1]]$layers[[1]]$aes_params$size <-
    ifelse(pointIsHighlighted, highlightPointSize, pointSize)

  # labels
  pUMAP[[1]]$layers[[2]]$aes_params$alpha <-
    ifelse(labelIsHighlighted, highlightLabelAlpha, labelAlpha)
  pUMAP[[1]]$layers[[2]]$aes_params$size <-
    ifelse(labelIsHighlighted, highlightLabelSize, labelSize)

  ### Add title
  pUMAP <-
    pUMAP +
    labs(title = 'Gene Set Clustering')

  if (interactive) {
    hoverInfo <- Seurat::FetchData(geneSetsList@cluster, vars = hoverVars)
    colnames(hoverInfo)[colnames(hoverInfo) == 'ident'] <- 'cluster'
    colnames(hoverInfo)[colnames(hoverInfo) == 'nCount_GSEA'] <- 'genes'
    return(Seurat::HoverLocator(plot = pUMAP[[1]], information = hoverInfo))
  }

  return(pUMAP)
}
