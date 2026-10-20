#' Plot enrichment as lollipop plot
#'
#' @description Creates a lollipop plot showing the summarized enrichment and size of each gene-set cluster.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param plotFrom Name of the column of the 'enrichmentScore' slot to visualize, e.g. 'nlog10pval' or 'enrichment_score'. Only numeric columns that are annotated (i.e. not all NA) can be plotted - which columns are available depends on what was supplied to 'geneSetsDF' in 'initializeList()'. Only one column is plotted; if several are supplied the first plottable one is used (the default therefore falls back to 'enrichment_score' when no p-values were supplied).
#' @param summaryStat How the values of the gene sets within a cluster are summarized into the value shown for each cluster. One of 'median' (default), 'max' or 'mean'. Note that 'max' is the largest signed value, so for signed scores (e.g. NES) a strongly negative cluster will be summarized by its least negative gene set.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the clusters (at the resolution) to highlight (in bold) in the plot.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plotEnrichments(exampleGsList)
#' plotEnrichments(exampleGsList, summaryStat = 'mean', highlightClusterNo = 0)
plotEnrichments <- function(
    geneSetsList,
    resolution = NULL,
    plotFrom = c('nlog10pval', 'enrichment_score'),
    summaryStat = c('median', 'max', 'mean'),
    removeClusterId = FALSE,
    highlightClusterNo = NULL
) {
  # Check input
  checkGsList(geneSetsList)

  summaryStat <- match.arg(summaryStat)

  plotFrom <- validatePlotFrom(
    geneSetsList = geneSetsList,
    plotFrom = plotFrom
  )

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
          size = genesets
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
        legend.position = "bottom",
        legend.box = "vertical" # side by side, the legends get clipped
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
