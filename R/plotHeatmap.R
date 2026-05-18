#' Plot enrichment as heatmap
#'
#' @description Creates a heatmap showing the mean enrichment of each gene-set cluster in each source (useful for datasets with many data sources).
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highlight in the UMAP
#' @param transpose A logical indicating whether to transpose the heatmap (switch rows and cols)
#' @param cutree_rows pheatmap control parameter. Number of clusters the rows are divided into
#' @param cutree_cols pheatmap control parameter. Number of clusters the columns are divided into
#' @param drawPlot Whether to show the plot. Default to TRUE.
#' @return A list of ggplot objects.
#' @export
#'
#' @examples plotHeatmap(exampleGsList)
plotHeatmap <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    highlightClusterNo = NULL,
    transpose = FALSE,
    cutree_rows = NA,
    cutree_cols = NA,
    drawPlot = TRUE
) {
  scaleCols = FALSE
  # Check input
  if (!inherits(geneSetsList, "gsList")) {
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
  if( all(is.na( geneSetsList@enrichmentScore$enrichment_score  ) ) ) {
    stop(paste0(
      'Cannot plot enrichment when enrichment scores are not annotated.\n',
      '  If you want to plot enrichment you need to redo the entire workflow and add the enrichment score.'
    ))
  }

  ### Extract clusters to highlight
  if(! is.null(highlightClusterNo)) {

    if( ! is.null(resolution)) {
      geneSetsList <- changeUsedResolution(
        geneSetsList = geneSetsList,
        resolution = resolution,
        removeClusterId = FALSE
      )
    }

    ### Extract cluster to highlight
    clusterPattern <- stringr::str_c('^',highlightClusterNo, ': ', collapse = '|')

    clusterNameToHighlight <-
      stringr::str_subset(
        string = levels(geneSetsList@cluster@active.ident),
        pattern = clusterPattern
      )

    if(removeClusterId) {
      clusterNameToHighlight <-
        stringr::str_remove(clusterNameToHighlight, "^[0-9]+: ")
    }
  } else {
    clusterNameToHighlight <- NULL
  }


  ### Update if not default
  if( ! is.null(resolution)) {
    geneSetsList <- changeUsedResolution(
      geneSetsList = geneSetsList,
      resolution = resolution,
      removeClusterId = removeClusterId
    )
  }

  ### Extract plot data
  plotData <-
    extractTable(
      geneSetsList = geneSetsList,
      removeClusterId = removeClusterId
    ) %>%
    dplyr::rename(
      mean = `enrichment mean`
    )

  ### Massage into matrix
  plotMatrix <-
    plotData %>%
    dplyr::select(source, cluster, mean) %>%
    tidyr::spread(key = cluster, value = mean) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('source') %>%
    as.matrix()

  # Clusters absent from a source produce NA after spread, which breaks
  # hclust downstream. Treat absent as zero enrichment.
  plotMatrix[is.na(plotMatrix)] <- 0

  if(transpose) {
    plotMatrix <- t(plotMatrix)
  }

  ### Determine colors
  if(any(!is.na(plotMatrix) & plotMatrix < 0)) {
    myColors <- colorRampPalette(c('blue','white',"red"))( 100 )
    myBreaks   <- seq(
      max(abs(plotMatrix), na.rm = T) * -1,
      max(abs(plotMatrix), na.rm = T),
      length.out = 99
    )
  } else {
    myColors <- colorRampPalette(c('white',"black"))( 100 )
    myBreaks   <- seq(
      min(plotMatrix, na.rm = T),
      max(plotMatrix, na.rm = T),
      length.out = 99
    )
  }

  ### Add legend title
  myLegendBreaks <- c(
    seq(
      max(abs(plotMatrix), na.rm = T) * -1,
      max(abs(plotMatrix), na.rm = T),
      length.out = 7
    )
  )
  myLegnedLabels <- round(myLegendBreaks, digits = 2)
  myLegnedLabels[length(myLegnedLabels)] <- paste0(
    'Enrich',
    '\n',
    myLegnedLabels[length(myLegnedLabels)],
    '\n'
  )

  ### Function to make rownames bold
  make_bold_names <- function(mat, rc_fun, rc_names) {
    bold_names <- rc_fun(mat)
    ids <- rc_names %>% match(rc_fun(mat))
    ids %>%
      purrr::walk(
        function(i)
          bold_names[i] <<-
          bquote(bold(.(rc_fun(mat)[i]))) %>%
          as.expression()
      )
    bold_names
  }


  ### Make plot
  tmp <- pheatmap::pheatmap(
    mat = plotMatrix,
    cluster_rows = T,
    cluster_cols = T,
    clustering_distance_rows = 'euclidean',
    clustering_distance_cols = 'euclidean',
    scale = ifelse(scaleCols, 'column','none'), # alternative 'none','column'
    kmeans_k = NA,     # When NA no aggregation is done
    show_rownames = T,
    show_colnames = T,
    color = myColors,
    breaks = myBreaks,
    legend_breaks = myLegendBreaks,
    legend_labels = myLegnedLabels,
    cutree_rows = cutree_rows, # set to NA to turn off
    cutree_cols = cutree_cols, # set to NA to turn off
    labels_row = make_bold_names(plotMatrix, rownames, clusterNameToHighlight),
    labels_col = make_bold_names(plotMatrix, colnames, clusterNameToHighlight),
    angle_col = 315,
    silent = ! drawPlot
  )
  if(!drawPlot) {
    return(tmp)
  } else {
    return(NULL)
  }
}
