#' Plot proportions of gene sets in each cluster for each source
#'
#' @description Graphs a bar plot showing the proportion of gene sets within
#' each cluster that comes from each of the two sources in the analyzed data.
#' Gene-sets found in both sources contribute to the fraction of both sources.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highligthClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highligth in the umap
#'
#' @return A list of ggplot objects.
#' @export
#'
#' @examples
#' plotProportions(exampleGsList)
plotProportions <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    highligthClusterNo = NULL
)  {
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
  if (!class(geneSetsList) == "gsList") {
    stop("'geneSetsList' should be a gsList object")
  }

  nCond <- countConditions(geneSetsList)
  # if( nCond != '2' ) {
  #   stop('plotProportions can only be used on data with two conditions')
  # }

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

  df_cluster <-
    geneSetsList@cluster@active.ident %>%
    tibble::enframe() %>%
    dplyr::rename("cluster" = "value", "pathway" = "name")

  df <- geneSetsList@enrichmentScore %>%
    dplyr::left_join(df_cluster, by = c("pathway"))

  if(removeClusterId) {
    df$cluster <- stringr::str_remove(df$cluster, "^[0-9]+: ")
  }

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

  df$cluster <- factor(
    df$cluster,
    levels = fctOrder$cluster
  )

  p <-
    df%>%
    dplyr::group_by(pathway) %>%
    dplyr::mutate(
      n_cond = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      source = dplyr::case_when(
        n_cond == 2 & nCond == 2 ~ 'Both',
        TRUE ~ source
      ),
      source = dplyr::case_when(
        n_cond >= 3 & nCond >= 3 ~ 'Multiple',
        TRUE ~ source
      )
    ) %>%
    ggplot(aes(x = cluster, fill = source)) +
    geom_bar(position = "fill") +
    labs(y = "Fraction of Gene Sets", x = "Gene Set Cluster", fill = "Source") +
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
  if(! is.null(highligthClusterNo)) {
    suppressWarnings(
      p <-
        p +
        theme(axis.text.y=element_text(
          face=boldXaxisLabels(
            src = df$cluster,
            clusterNameToHighligth
          )
        ))
    )

  }


  return(p)
}
