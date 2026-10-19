#' Extract summary table of clustering
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution Numeric. Resolution to extract the table for. Default
#'   `NULL` uses the active resolution: the optimal resolution found by
#'   `clusterGeneSets()`, or the one last set with `updateResolution()` or
#'   `updateClusterLabels()`.
#' @param summarize Logical indicating whether to summarize the table by
#'   cluster (and source). Default is TRUE.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed
#'   from each cluster label. Default: TRUE.
#'
#' @return A data.frame. With `summarize = TRUE`, one row per cluster and
#'   source with the mean and standard deviation of the enrichment score
#'   (ignoring missing values) and the number of gene sets. With
#'   `summarize = FALSE`, one row per gene set and source with all columns of
#'   the 'enrichmentScore' slot plus its `cluster`; if a cluster label was
#'   changed with `updateClusterLabels()`, the original label is given in
#'   `initial_cluster_label`.
#' @export
#'
#' @examples
#' extractTable(exampleGsList)
#' head(extractTable(exampleGsList, resolution = 0.6, summarize = FALSE))
extractTable <- function(
    geneSetsList,
    resolution = NULL,
    summarize = TRUE,
    removeClusterId = TRUE
) {

  checkGsList(geneSetsList)

  resolutionColumn <- resolveResolutionColumn(geneSetsList, resolution)

  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = resolution,
    removeClusterId = removeClusterId
  )

  ### Get cluster labels
  df_cluster <- data.frame(
    pathway = names(geneSetsList@cluster@active.ident),
    cluster = as.character(geneSetsList@cluster@active.ident),
    stringsAsFactors = FALSE
  )

  ### Add original label
  labelsDefault <- geneSetsList@labels[[resolutionColumn]]$labelsDefault
  if (!summarize && 'org_label' %in% colnames(labelsDefault)) {
    clusterIds <- as.character(
      geneSetsList@cluster@meta.data[[resolutionColumn]]
    )
    names(clusterIds) <- rownames(geneSetsList@cluster@meta.data)

    originalLabels <- labelsDefault %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(weighted_score, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()

    df_cluster$initial_cluster_label <- originalLabels$org_label[match(
      clusterIds[df_cluster$pathway],
      as.character(originalLabels$cluster)
    )]
  }

  ### Keep the cluster order of the plots (numeric cluster id)
  clusterOrder <- levels(geneSetsList@cluster@active.ident)

  ### Extract info
  clusterScores <-
    geneSetsList@enrichmentScore %>%
    dplyr::inner_join(
      y = df_cluster,
      by = "pathway"
    ) %>%
    dplyr::arrange(match(cluster, clusterOrder))

  ### Summarize
  if( summarize ) {
    clusterScores <-
      clusterScores %>%
      dplyr::group_by(source, cluster) %>%
      dplyr::summarise(
        `enrichment mean` = round(
          summarizeScores(enrichment_score, mean), digits = 2
        ),
        `enrichment sd` = round(
          summarizeScores(enrichment_score, stats::sd), digits = 2
        ),
        genesets = dplyr::n(),
        .groups = 'drop'
      ) %>%
      dplyr::arrange(match(cluster, clusterOrder))
  }

  return(as.data.frame(clusterScores))
}
