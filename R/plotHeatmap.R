#' Plot enrichment as heatmap
#'
#' @description Creates a heatmap summarizing the enrichment of each gene-set
#'   cluster in each source (useful for datasets with many data sources).
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param plotFrom Name of the column of the 'enrichmentScore' slot to visualize, e.g. 'nlog10pval' or 'enrichment_score'. Only numeric columns that are annotated (i.e. not all NA) can be plotted - which columns are available depends on what was supplied to 'geneSetsDF' in 'initializeList()'. Only one column is plotted; if several are supplied the first plottable one is used (the default therefore falls back to 'enrichment_score' when no p-values were supplied).
#' @param summaryStat How the values of the gene sets within a cluster are summarized into the value shown in each heatmap cell. One of 'median' (default), 'max' or 'mean'. Note that 'max' is the largest signed value, so for signed scores (e.g. NES) a strongly negative cluster will be summarized by its least negative gene set.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the clusters (at the resolution) to highlight (in bold) in the plot.
#' @param transpose A logical indicating whether to transpose the heatmap (switch rows and cols)
#' @param cutree_rows pheatmap control parameter. Number of clusters the rows are divided into
#' @param cutree_cols pheatmap control parameter. Number of clusters the columns are divided into
#' @param drawPlot Whether to draw the plot. Default to TRUE.
#' @param scaleCols Whether to z-score each cluster across the sources before
#'   plotting. One of `"auto"` (default), `TRUE` or `FALSE`. `"auto"` scales
#'   only when there are at least three sources: with one source there is
#'   nothing to compare, and with two sources the z-scores can only be
#'   -0.71, 0 or 0.71. Scaling is always across sources, also when
#'   `transpose = TRUE`.
#' @return The pheatmap object, invisibly. Use `drawPlot = FALSE` to get it
#'   without drawing, e.g. to arrange it with other plots.
#'
#' @details Cells of clusters that have no gene set from a source are shown in
#'   gray. Rows and columns are ordered by hierarchical clustering when there
#'   are at least two sources (absent cells count as 0 for the clustering).
#' @export
#'
#' @examples
#' plotHeatmap(exampleGsList)
#' plotHeatmap(exampleGsList, plotFrom = 'enrichment_score', summaryStat = 'mean')
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
    scaleCols = "auto"
) {
  # Check input
  checkGsList(geneSetsList)

  summaryStat <- match.arg(summaryStat)

  plotFrom <- validatePlotFrom(
    geneSetsList = geneSetsList,
    plotFrom = plotFrom
  )

  if (identical(scaleCols, "auto")) {
    scaleCols <- countConditions(geneSetsList) == '3+'
  }
  if (!is.logical(scaleCols) || length(scaleCols) != 1L || is.na(scaleCols)) {
    stop("'scaleCols' should be \"auto\", TRUE or FALSE.")
  }

  clusterNameToHighlight <- highlightedClusterNames(
    geneSetsList = geneSetsList,
    resolution = resolution,
    highlightClusterNo = highlightClusterNo,
    removeClusterId = removeClusterId
  )

  ### Extract plot data
  plotData <-
    summarizeEnrichment(
      geneSetsList = geneSetsList,
      resolution = resolution,
      scoreColumn = plotFrom,
      summaryStat = summaryStat,
      removeClusterId = removeClusterId
    )

  # Data without a 'source' column has a single, unnamed source
  plotData$source[is.na(plotData$source)] <- 'All gene sets'

  ### Massage into matrix (sources x clusters)
  clusterOrder <- unique(plotData$cluster)
  plotMatrix <-
    plotData %>%
    dplyr::select(source, cluster, enrichmentSummary) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = enrichmentSummary) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('source') %>%
    as.matrix()
  plotMatrix <- plotMatrix[, clusterOrder, drop = FALSE]

  ### Normalize each cluster across sources
  # Done here rather than through pheatmap's own 'scale' argument because the
  # colors, breaks and legend below are all derived from 'plotMatrix' -- if
  # pheatmap scaled internally those would stay on the raw scale.
  if(scaleCols) {
    isMissing <- is.na(plotMatrix)
    plotMatrix <- scale(plotMatrix, center = TRUE, scale = TRUE)

    # Clusters without any variation across sources become NaN - they are
    # flat, so 0 is the sensible z-score for every one of their cells
    plotMatrix[is.nan(plotMatrix) & !isMissing] <- 0

    attr(plotMatrix, 'scaled:center') <- NULL
    attr(plotMatrix, 'scaled:scale')  <- NULL
  }

  if(transpose) {
    plotMatrix <- t(plotMatrix)
  }

  ### Hierarchical clustering of rows and columns
  # Done on a copy where absent cells are 0, so the heatmap itself can show
  # them as missing.
  clusteringMatrix <- plotMatrix
  clusteringMatrix[is.na(clusteringMatrix)] <- 0
  doHCluster <- countConditions(geneSetsList) != '1'
  clusterRows <- FALSE
  clusterCols <- FALSE
  if (doHCluster && nrow(plotMatrix) > 1) {
    clusterRows <- stats::hclust(stats::dist(clusteringMatrix))
  }
  if (doHCluster && ncol(plotMatrix) > 1) {
    clusterCols <- stats::hclust(stats::dist(t(clusteringMatrix)))
  }

  ### Determine colors
  if(any(!is.na(plotMatrix) & plotMatrix < 0)) {
    # Signed values (e.g. NES, or anything z-scored) get a diverging palette
    # centered on zero
    myColors <- colorRampPalette(c('royalblue','white',"#e34632"))( 100 )
    breakRange <- c(-1, 1) * max(abs(plotMatrix), na.rm = TRUE)
  } else {
    # Non-negative values (e.g. unscaled nlog10pval) get a sequential palette
    myColors <- colorRampPalette(c('white',"black"))( 100 )
    breakRange <- range(plotMatrix, na.rm = TRUE)
  }

  # A constant matrix gives a zero-width range; widen it so the breaks are
  # unique
  if (diff(breakRange) == 0) {
    breakRange <- breakRange + c(-1, 1) * max(abs(breakRange[1]), 1)
  }

  myBreaks <- seq(breakRange[1], breakRange[2], length.out = 101)

  ### Add legend title
  # Taken from the same range as 'myBreaks' - otherwise a matrix without
  # negative values gets a legend spanning values its colors never cover
  myLegendBreaks <- seq(breakRange[1], breakRange[2], length.out = 7)
  myLegendLabels <- round(myLegendBreaks, digits = 2)
  myLegendLabels[length(myLegendLabels)] <- paste0(
    ifelse(scaleCols, 'Z-score', 'Raw'),
    '\n',
    myLegendLabels[length(myLegendLabels)],
    '\n'
  )

  ### Title naming the value shown
  # Uses the column 'validatePlotFrom()' settled on, so it also tells the user
  # when the default fell back from 'nlog10pval' to 'enrichment_score'. The
  # scaling is left to the legend title, which already says 'Z-score'/'Raw'.
  plotTitle <- enrichmentSummaryLabel(summaryStat, plotFrom)

  ### Make plot
  tmp <- pheatmap::pheatmap(
    mat = plotMatrix,
    cluster_rows = clusterRows,
    cluster_cols = clusterCols,
    scale = 'none', # 'scaleCols' is already applied to 'plotMatrix' above
    kmeans_k = NA,     # When NA no aggregation is done
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = myColors,
    breaks = myBreaks,
    na_col = 'grey85',
    legend_breaks = myLegendBreaks,
    legend_labels = myLegendLabels,
    cutree_rows = cutree_rows, # set to NA to turn off
    cutree_cols = cutree_cols, # set to NA to turn off
    labels_row = boldNames(rownames(plotMatrix), clusterNameToHighlight),
    labels_col = boldNames(colnames(plotMatrix), clusterNameToHighlight),
    angle_col = 315,
    main = plotTitle,
    silent = ! drawPlot
  )

  return(invisible(tmp))
}

#' Plotmath labels with selected names in bold
#'
#' @param names Character vector of row or column names.
#' @param toBold Names to show in bold.
#' @return An expression vector usable as pheatmap labels.
#' @keywords internal
boldNames <- function(names, toBold) {
  labels <- lapply(names, function(name) {
    if (name %in% toBold) {
      bquote(bold(.(name)))
    } else {
      bquote(.(name))
    }
  })
  as.expression(labels)
}
