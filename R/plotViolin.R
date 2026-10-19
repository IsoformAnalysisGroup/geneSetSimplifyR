#' Violin plot of clustering
#'
#' @description Creates a violin plot of the values of the gene sets in each
#'   cluster, by default their -log10 p-values (or enrichment scores when no
#'   p-values were supplied). Clusters are ordered by their median value.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param plotFrom Name of the column of the 'enrichmentScore' slot to visualize, e.g. 'nlog10pval' or 'enrichment_score'. Only numeric columns that are annotated (i.e. not all NA) can be plotted - which columns are available depends on what was supplied to 'geneSetsDF' in 'initializeList()'. Only one column is plotted; if several are supplied the first plottable one is used (the default therefore falls back to 'enrichment_score' when no p-values were supplied).
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the clusters (at the resolution) to highlight (in bold) in the plot.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plotViolin(exampleGsList, highlightClusterNo = c(0, 1))
#' plotViolin(exampleGsList, plotFrom = 'enrichment_score', removeClusterId = TRUE)
plotViolin <- function(
    geneSetsList,
    resolution = NULL,
    plotFrom = c('nlog10pval', 'enrichment_score'),
    removeClusterId = FALSE,
    highlightClusterNo = NULL
) {
  # Check input
  checkGsList(geneSetsList)

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

  ### Set the cluster resolution and labels to plot
  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = resolution,
    removeClusterId = removeClusterId
  )

  ### One row per gene set and source, with its cluster label
  clusterIdents <- geneSetsList@cluster@active.ident
  df <- geneSetsList@enrichmentScore
  df$clusterLabel <- clusterIdents[match(df$pathway, names(clusterIdents))]
  df$value <- df[[plotFrom]]
  df <- df[!is.na(df$value) & !is.na(df$clusterLabel), , drop = FALSE]

  ### Order the clusters by their median value
  df$clusterLabel <- droplevels(df$clusterLabel)
  df$clusterLabel <- stats::reorder(df$clusterLabel, df$value, FUN = median)

  valueLabel <- sub("^ ", "", enrichmentSummaryLabel('', plotFrom))

  # A violin needs at least two values; single gene sets are still drawn as
  # points
  violinData <- df %>%
    dplyr::group_by(clusterLabel, source) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  nCond <- countConditions(geneSetsList)

  if ( nCond != '1' ) {
    p1 <-
      df %>%
      ggplot(aes(
        x = clusterLabel,
        y = value,
        fill = source
      )) +
      violinWithMedian(violinData) +
      geom_point(
        alpha = 1/3,
        size = 0.75,
        position = position_jitterdodge(
          jitter.width = 0.2,
          jitter.height = 0,
          dodge.width = 0.9,
          seed = 1
        )
      ) +
      labs(x = "Gene Set Cluster",
           y = valueLabel,
           fill = "Gene set source") +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      coord_flip()
  }
  if ( nCond == '1' ) {
    df <- df %>%
      dplyr::group_by(clusterLabel) %>%
      dplyr::mutate(median_value = median(value)) %>%
      dplyr::ungroup()
    violinData <- df %>%
      dplyr::group_by(clusterLabel) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::ungroup()

    p1 <-
      df %>%
      ggplot(aes(
        x = clusterLabel,
        y = value,
        fill = median_value
      )) +
      violinWithMedian(violinData) +
      geom_point(
        alpha = 1/3,
        size = 0.75,
        position = position_jitter(width = 0.2, height = 0, seed = 1)
      ) +
      labs(x = "Gene Set Cluster",
           y = valueLabel,
           fill = paste("Median", valueLabel)) +
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
            src = df$clusterLabel,
            clusterNameToHighlight
          )
        ))
    )
  }

  return(p1)
}

#' Violin layer with a line at the median
#'
#' ggplot2 4.0 replaced `draw_quantiles` by `quantiles`/`quantile.linetype`;
#' use whichever the installed version supports.
#' @param data Data for the layer.
#' @return A ggplot2 layer.
#' @keywords internal
violinWithMedian <- function(data) {
  if (utils::packageVersion("ggplot2") >= "4.0.0") {
    geom_violin(data = data, scale = 'width',
                quantiles = 0.5, quantile.linetype = 1)
  } else {
    geom_violin(data = data, scale = 'width', draw_quantiles = 0.5)
  }
}
