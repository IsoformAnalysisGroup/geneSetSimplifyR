#' Plot clusters as a bar chart
#'
#' @description Creates a bar chart showing the size of each gene-set cluster colored by condition.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highligthClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highligth in the umap
#'
#' @return A list of ggplot objects.
#' @export
#'
#' @examples plotBarchart(geneSetsList)
plotBarchart <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    highligthClusterNo = NULL
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

  ### Extract clusters to highligth
  if(! is.null(highligthClusterNo)) {

    if( ! is.null(resolution)) {
      geneSetsList <- changeUsedResolution(
        geneSetsList = geneSetsList,
        resolution = resolution,
        removeClusterId = FALSE
      )
    }

    ### Extact cluster to highligth
    clusterPattern <- stringr::str_c('^',highligthClusterNo, ': ', collapse = '|')

    clusterNameToHighligth <-
      stringr::str_subset(
        string = levels(geneSetsList@cluster@active.ident),
        pattern = clusterPattern
      )

    if(removeClusterId) {
      clusterNameToHighligth <-
        stringr::str_remove(clusterNameToHighligth, "^[0-9]+: ")
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
    extractTable(
      geneSetsList,
      removeClusterId = removeClusterId,
      summarize = T
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
        y = "Gene Set",
        color = "Gene set source",
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
        y = "Gene Set",
        color = "Gene set source",
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      coord_flip()
  }

  ### Make highlighted clusters bold
  if(! is.null(highligthClusterNo)) {
    suppressWarnings(
      p1 <-
        p1 +
        theme(axis.text.y=element_text(
          face=boldXaxisLabels(
            src = plotData$cluster,
            clusterNameToHighligth
          )
        ))
    )

  }

  return(p1)
}
