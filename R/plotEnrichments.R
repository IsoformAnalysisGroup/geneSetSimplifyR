#' Plot enrichment as lollipop plot
#'
#' @description Creates a lollipop plot showing the summarized enrichment and size of each gene-set cluster.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param plotFrom Name of the column of the 'enrichmentScore' slot to visualize, e.g. 'enrichment_score' or 'nlog10pval'. Only numeric columns that are annotated (i.e. not all NA) can be plotted - which columns are available depends on what was supplied to 'geneSetsDF' in 'initializeList()'. Only one column is plotted; if several are supplied the first plottable one is used.
#' @param summaryStat How the values of the gene sets within a cluster are summarized into the value shown for each cluster. One of 'mean' (default), 'median' or 'max'. Note that 'max' is the largest signed value, so for signed scores (e.g. NES) a strongly negative cluster will be summarized by its least negative gene set.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highlight in the UMAP
#'
#' @return A ggplot object.
#' @export
#'
#' @examples plotEnrichments(exampleGsList)
plotEnrichments <- function(
    geneSetsList,
    resolution = NULL,
    plotFrom = c('enrichment_score', 'nlog10pval'),
    summaryStat = c('mean', 'median', 'max'),
    removeClusterId = FALSE,
    highlightClusterNo = NULL
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

  yAxisLabel <- enrichmentSummaryLabel(summaryStat, plotFrom)

  plotData$cluster <- factor(
    plotData$cluster,
    levels = rev(unique(plotData$cluster))
  )


  nCond <- countConditions(geneSetsList)

  if ( nCond != '1' ) {
    p1 <-
      plotData %>%
      ggplot(aes(x=cluster, y=enrichmentSummary, group=source)) +
      geom_linerange(
        mapping = aes(
          x=cluster,
          ymin=0,
          ymax=enrichmentSummary,
          color=source
        ),
        position=position_dodge(width = 0.6),
        show.legend = FALSE # the dots already carry the color legend
      ) +
      geom_point(
        mapping = aes(
          color=source,
          size = genesets,
        ),
        position=position_dodge(width = 0.6)
      ) +
      scale_radius(
        limits = range(plotData$genesets),
        breaks = round(seq(
          min(plotData$genesets),
          max(plotData$genesets),
          length.out = 4
        ))
      ) +
      labs(
        x = "Gene Set Cluster",
        y = yAxisLabel,
        color = "Gene set source",
        size = "Gene Sets"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
      coord_flip()
  }

  if ( nCond == '1' ) {
    p1 <-
      plotData %>%
      ggplot(aes(x=cluster, y=enrichmentSummary, group=source)) +
      geom_linerange(
        mapping = aes(
          x=cluster,
          ymin=0,
          ymax=enrichmentSummary
        ),
        color="#595959",
        position=position_dodge(width = 0.6)
      ) +
      geom_point(
        mapping = aes(
          size = genesets
        ),
        color="#595959",
        position=position_dodge(width = 0.6)
      ) +
      scale_radius(
        limits = range(plotData$genesets),
        breaks = round(seq(
          min(plotData$genesets),
          max(plotData$genesets),
          length.out = 4
        ))
      ) +
      labs(
        x = "Gene Set Cluster",
        y = yAxisLabel,
        size = "Gene Sets"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      coord_flip()
  }

  ### Make highlighted clusters bold
  if(! is.null(highlightClusterNo)) {
    suppressWarnings(
      p1 <-
        p1 +
        theme(axis.text.y=element_text(
          face=boldXaxisLabels(
            src = plotData$cluster,
            clusterNameToHighlight
          )
        ))
    )

  }

  return(p1)
}
