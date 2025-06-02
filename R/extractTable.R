#' Extract summary table of clustering
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution Numeric. Resolution to plot the table for.
#' @param summarize Logic indicating whether to summarize table by cluster (and source). Default is TRUE.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#'
#' @return Table with clustering summary
#' @export
#'
#' @examples
#' extractTable(exampleGsList, 1.1)
extractTable <- function(
    geneSetsList,
    resolution = NULL,
    summarize = TRUE,
    removeClusterId = TRUE
) {

  # Check input
  if (!class(geneSetsList) == "gsList") {
    stop("'geneSetsList' should be a gsList object")
  }

  if (length(geneSetsList@cluster) == 0) {
    stop("No clustering performed for 'geneSetsList'.
         Please run 'clusterGeneSets()' and 'labelClusters()' before plotting.")
  }

  if (length(geneSetsList@labels) == 0) {
    stop("No labels created for 'geneSetsList'.
         Please run 'labelClusters()' before plotting.")
  }

  ### Update if not default
  if( ! is.null(resolution)) {
    geneSetsList <- changeUsedResolution(
      geneSetsList = geneSetsList,
      resolution = resolution,
      removeClusterId = FALSE
    )
  }

  if( ! is.null(resolution)) {
    resolutionColumn <-
      paste0("GSEA_snn_res.", resolution)
  } else {
    resolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]
  }

  ### Get cluster labels
  df_cluster <-
    geneSetsList@cluster@active.ident %>%
    tibble::enframe() %>%
    dplyr::rename("cluster" = "value", "pathway" = "name")

  ### Add original label
  t1 <- ! summarize
  t2 <- 'org_label' %in% colnames(geneSetsList@labels[[resolutionColumn]])
  if( t1 & t2 ) {
    df_label <-
      geneSetsList@labels[[resolutionColumn]]$labelsDefault %>%
      dplyr::mutate(
        cluster = stringr::str_c(
          cluster,
          ': ',
          label
        )
      ) %>%
      dplyr::select(
        cluster,
        org_label
      )

    df_cluster$initial_cluster_label <-
      df_label$org_label[match(
        df_cluster$cluster, df_label$cluster
      )]
  }

  ### Extract info
  clusterScores <-
    geneSetsList@enrichmentScore %>%
    dplyr::inner_join(
      x = .,
      y = df_cluster,
      by = "pathway"
    ) %>%
    dplyr::arrange(cluster)

  ### Summarize
  if( summarize ) {
    clusterScores <-
      clusterScores %>%
      dplyr::group_by(source, cluster) %>%
      dplyr::summarise(
        `enrichment mean` = round(mean(enrichment_score), digits = 2),
        `enrichment sd` = round(sd(enrichment_score), digits = 2),
        genesets = dplyr::n(),
        .groups = 'drop'
      ) %>%
      dplyr::arrange(cluster)
  }

  if (removeClusterId) {
    clusterScores$cluster <- stringr::str_remove(clusterScores$cluster, "^[0-9]+: ")
  }

  return(clusterScores)
}
