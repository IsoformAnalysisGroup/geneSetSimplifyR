#' Plot clusters as a bar chart
#'
#' @description Creates a bar chart showing the number of gene sets in each
#'   cluster, per source when there are several.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the clusters (at the resolution) to highlight (in bold) in the plot.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plotBarchart(exampleGsList)
#' plotBarchart(exampleGsList, resolution = 0.6, highlightClusterNo = 1)
plotBarchart <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    highlightClusterNo = NULL
) {
  # Check input
  checkGsList(geneSetsList)

  clusterNameToHighlight <- highlightedClusterNames(
    geneSetsList = geneSetsList,
    resolution = resolution,
    highlightClusterNo = highlightClusterNo,
    removeClusterId = removeClusterId
  )

  ### Extract plot data
  plotData <-
    extractTable(
      geneSetsList,
      resolution = resolution,
      removeClusterId = removeClusterId,
      summarize = TRUE
    )

  plotData$cluster <- factor(
    plotData$cluster,
    levels = rev(unique(plotData$cluster))
  )

  nCond <- countConditions(geneSetsList)

  if ( nCond != '1' ) {
    p1 <-
      plotData %>%
      ggplot(aes(x=cluster, y=genesets, fill=source)) +
      geom_bar(stat="identity", position=position_dodge2(preserve = "single")) +
      labs(
        x = "Gene Set Cluster",
        y = "Number of gene sets",
        fill = "Gene set source"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      coord_flip()
  }

  if ( nCond == '1' ) {
    p1 <-
      plotData %>%
      ggplot(aes(x=cluster, y=genesets)) +
      geom_bar(stat="identity", position=position_dodge2(preserve = "single")) +
      labs(
        x = "Gene Set Cluster",
        y = "Number of gene sets"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      coord_flip()
  }

  ### Make highlighted clusters bold
  if(length(clusterNameToHighlight)) {
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
