#' Plot proportions of gene sets in each cluster for each source
#'
#' @description Graphs a bar plot showing the proportion of gene sets within
#' each cluster that come from each source in the analyzed data. A gene set
#' found in several sources is counted once, as 'Both' (two sources) or
#' 'Multiple' (three or more). Needs at least two sources.
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
#' plotProportions(exampleGsList)
plotProportions <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    highlightClusterNo = NULL
)  {
  # Check input
  checkGsList(geneSetsList)

  nCond <- countConditions(geneSetsList)
  if( nCond == '1' ) {
    stop(
      "plotProportions() compares sources and needs gene sets from at least ",
      "two sources (the 'source' column of 'geneSetsDF' in 'initializeList()')."
    )
  }

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

  clusterIdents <- geneSetsList@cluster@active.ident
  df <- geneSetsList@enrichmentScore
  df$cluster <- as.character(clusterIdents[match(df$pathway, names(clusterIdents))])

  groupSize <-
    df %>%
    dplyr::group_by(source, cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      .groups = 'drop_last'
    ) %>%
    dplyr::summarise(
      frac_cluster = sum(n > 0) / length(unique(df$cluster)),
      mean_size = mean(n),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(dplyr::desc(frac_cluster), dplyr::desc(mean_size))

  fctOrder <-
    df %>%
    dplyr::group_by(cluster, source) %>%
    dplyr::summarise(
      n = dplyr::n(),
      .groups = 'drop_last'
    ) %>%
    dplyr::mutate(
      frac = n / sum(n)
    ) %>%
    dplyr::filter(
      source == groupSize$source[1]
    ) %>%
    dplyr::arrange(dplyr::desc(frac))

  allClusters <- unique(as.character(df$cluster))
  orderedLevels <- c(
    as.character(fctOrder$cluster),
    setdiff(allClusters, as.character(fctOrder$cluster))
  )
  df$cluster <- factor(
    df$cluster,
    levels = orderedLevels
  )

  ### Gene sets found in several sources are counted once, as 'Both' (two
  ### sources) or 'Multiple' (three or more sources)
  plotData <-
    df %>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(
      n_cond = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      source = dplyr::case_when(
        n_cond >= 2 & nCond == '2' ~ 'Both',
        n_cond >= 2 ~ 'Multiple',
        TRUE ~ source
      )
    ) %>%
    dplyr::distinct(pathway, .keep_all = TRUE)

  p <-
    plotData %>%
    ggplot(aes(x = cluster, fill = source)) +
    geom_bar(position = "fill") +
    labs(y = "Fraction of Gene Sets", x = "Gene Set Cluster", fill = "Gene set source") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_flip() +
    theme(
      legend.position = "bottom"
    )

  ### Make highlighted clusters bold
  if(length(clusterNameToHighlight)) {
    suppressWarnings(
      p <-
        p +
        theme(axis.text.y=element_text(
          face=boldXaxisLabels(
            src = df$cluster,
            clusterNameToHighlight
          )
        ))
    )

  }


  return(p)
}
