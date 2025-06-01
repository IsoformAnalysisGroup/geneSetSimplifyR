#' Plot enrichment as lollipop plot
#'
#' @description Creates a lollipop plot of showing mean enrichment and size of each gene-set cluster.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highligthClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highligth in the umap
#'
#' @return A list of ggplot objects.
#' @export
#'
#' @examples plotEnrichments(geneSetsList)
plotEnrichments <- function(
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
  if( all(is.na( geneSetsList@enrichmentScore$enrichment_score  ) ) ) {
    stop(paste0(
      'Cannot plot enrichment when enrichment scores are not annoated.\n',
      '  If you want to plot enrichment you need to redo the entire workflow and add the enrichment score.'
    ))
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
    extractTable(geneSetsList, removeClusterId = removeClusterId) %>%
    dplyr::rename(
      mean = `enrichment mean`
    )

  plotData$cluster <- factor(
    plotData$cluster,
    levels = rev(unique(plotData$cluster))
  )


  nCond <- countConditions(geneSetsList)

  if ( nCond != '1' ) {
    p1 <-
      plotData %>%
      ggplot(aes(x=cluster, y=mean, group=source)) +
      geom_linerange(
        mapping = aes(
          x=cluster,
          ymin=0,
          ymax=mean
        ),
        color="#595959",
        position=position_dodge(width = 0.9)
      ) +
      geom_point(
        mapping = aes(
          color=source,
          size = genesets,
        ),
        position=position_dodge(width = 0.9)
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
        y = "Mean Enrichment",
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
      ggplot(aes(x=cluster, y=mean, group=source)) +
      geom_linerange(
        mapping = aes(
          x=cluster,
          ymin=0,
          ymax=mean
        ),
        color="#595959",
        position=position_dodge(width = 0.9)
      ) +
      geom_point(
        mapping = aes(
          size = genesets
        ),
        color="#595959",
        position=position_dodge(width = 0.9)
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
        y = "Mean Enrichment",
        size = "Gene Sets"
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
