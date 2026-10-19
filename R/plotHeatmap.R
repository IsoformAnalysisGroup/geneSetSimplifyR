#' Plot enrichment as heatmap
#'
#' @description Creates a heatmap summarizing the enrichment of each gene-set cluster in each source (useful for datasets with many data sources).
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param plotFrom Name of the column of the 'enrichmentScore' slot to visualize, e.g. 'nlog10pval' or 'enrichment_score'. Only numeric columns that are annotated (i.e. not all NA) can be plotted - which columns are available depends on what was supplied to 'geneSetsDF' in 'initializeList()'. Only one column is plotted; if several are supplied the first plottable one is used (the default therefore falls back to 'enrichment_score' when no p-values were supplied).
#' @param summaryStat How the values of the gene sets within a cluster are summarized into the value shown in each heatmap cell. One of 'median' (default), 'max' or 'mean'. Note that 'max' is the largest signed value, so for signed scores (e.g. NES) a strongly negative cluster will be summarized by its least negative gene set.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highlight in the UMAP
#' @param transpose A logical indicating whether to transpose the heatmap (switch rows and cols)
#' @param cutree_rows pheatmap control parameter. Number of clusters the rows are divided into
#' @param cutree_cols pheatmap control parameter. Number of clusters the columns are divided into
#' @param drawPlot Whether to show the plot. Default to TRUE.
#' @param scaleCols Whether to perform standardization.
#' @return The pheatmap object if 'drawPlot' is FALSE, else NULL (the heatmap is drawn).
#' @export
#'
#' @examples plotHeatmap(exampleGsList)
plotHeatmap <- function(
    geneSetsList,
    resolution = NULL,
    plotFrom = c('nlog10pval', 'enrichment_score'),
    summaryStat = c('median', 'max', 'mean'),
    removeClusterId = FALSE,
    highlightClusterNo = NULL,
    transpose = FALSE,
    cutree_rows = NA,
    cutree_cols = NA,
    drawPlot = TRUE,
    scaleCols = TRUE
) {
  # Check input
  if (!inherits(geneSetsList, "gsList")) {
    stop("'geneSetsList' should be a gsList object")
  }

  summaryStat <- match.arg(summaryStat)

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
  plotFrom <- validatePlotFrom(
    geneSetsList = geneSetsList,
    plotFrom = plotFrom
  )

  if (geneSetsList@metadata$n_sources == 1) {
    doHCluster <- FALSE
  } else{
    doHCluster <- TRUE
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
    summarizeEnrichment(
      geneSetsList = geneSetsList,
      scoreColumn = plotFrom,
      summaryStat = summaryStat,
      removeClusterId = removeClusterId
    )

  ### Massage into matrix
  plotMatrix <-
    plotData %>%
    dplyr::select(source, cluster, enrichmentSummary) %>%
    tidyr::spread(key = cluster, value = enrichmentSummary) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('source') %>%
    as.matrix()

  # Clusters absent from a source produce NA after spread, which breaks
  # hclust downstream. Treat absent as zero enrichment.
  plotMatrix[is.na(plotMatrix)] <- 0

  if(transpose) {
    plotMatrix <- t(plotMatrix)
  }

  ### Normalize columns
  # Done here rather than through pheatmap's own 'scale' argument because the
  # colors, breaks and legend below are all derived from 'plotMatrix' -- if
  # pheatmap scaled internally those would stay on the raw enrichment scale.
  if(scaleCols) {
    plotMatrix <- scale(plotMatrix, center = TRUE, scale = TRUE)

    # Columns without any variation become NaN - they are flat, so 0 is the
    # sensible z-score for every one of their cells
    plotMatrix[is.nan(plotMatrix)] <- 0

    attr(plotMatrix, 'scaled:center') <- NULL
    attr(plotMatrix, 'scaled:scale')  <- NULL
  }

  ### Determine colors
  if(any(!is.na(plotMatrix) & plotMatrix < 0)) {
    # Signed values (e.g. NES, or anything z-scored) get a diverging palette
    # centered on zero
    myColors <- colorRampPalette(c('royalblue','white',"#e34632"))( 100 )
    breakRange <- c(-1, 1) * max(abs(plotMatrix), na.rm = T)
  } else {
    # Non-negative values (e.g. unscaled nlog10pval) get a sequential palette
    myColors <- colorRampPalette(c('white',"black"))( 100 )
    breakRange <- range(plotMatrix, na.rm = T)
  }

  myBreaks <- seq(breakRange[1], breakRange[2], length.out = 99)

  ### Add legend title
  # Taken from the same range as 'myBreaks' - otherwise a matrix without
  # negative values gets a legend spanning values its colors never cover
  myLegendBreaks <- seq(breakRange[1], breakRange[2], length.out = 7)
  myLegnedLabels <- round(myLegendBreaks, digits = 2)
  myLegnedLabels[length(myLegnedLabels)] <- paste0(
    ifelse(scaleCols, 'Z-score', 'Raw'),
    '\n',
    myLegnedLabels[length(myLegnedLabels)],
    '\n'
  )

  ### Title naming the value shown
  # Uses the column 'validatePlotFrom()' settled on, so it also tells the user
  # when the default fell back from 'nlog10pval' to 'enrichment_score'. The
  # scaling is left to the legend title, which already says 'Z-score'/'Raw'.
  plotTitle <- enrichmentSummaryLabel(summaryStat, plotFrom)

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
    cluster_rows = doHCluster,
    cluster_cols = doHCluster,
    clustering_distance_rows = 'euclidean',
    clustering_distance_cols = 'euclidean',
    scale = 'none', # 'scaleCols' is already applied to 'plotMatrix' above
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
    main = plotTitle,
    silent = ! drawPlot
  )
  if(!drawPlot) {
    return(tmp)
  } else {
    return(NULL)
  }
}
