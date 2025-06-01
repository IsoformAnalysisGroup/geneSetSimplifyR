#' Violin plot of clustering
#'
#' @description Creates a violin plot of the enrichment scores for each cluster.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highligthClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highligth in the umap
#'
#' @return A list of ggplot objects.
#' @importFrom dplyr mutate
#' @export
#'
#' @examples plotViolin(geneSetsList)
plotViolin <- function(
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

  ### Assign conditions
  nCond <- countConditions(geneSetsList)

  if ( nCond != '1' ) {

    df <-
      geneSetsList@cluster %>%
      #tidyseurat::mutate(label = geneSetsList@cluster@active.ident) %>%
      mutate(label = geneSetsList@cluster@active.ident) %>%
      dplyr::right_join(.,
                        geneSetsList@enrichmentScore,
                        by = c("gene_set" = "pathway")) %>%
      dplyr::group_by(label) %>%
      dplyr::mutate(median_enrichment_score = median(enrichment_score)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        id = stringr::str_c(source, label),
        label = reorder(label, enrichment_score, FUN = median)
      )

    if(removeClusterId) {
      df$label <- stringr::str_remove(df$label, "^[0-9]+: ")
    }


    p1 <-
      df %>%
      ggplot(aes(
        x = label,
        y = enrichment_score,
        fill = source
      )) +
      geom_violin(scale = 'width', draw_quantiles = 0.5) +
      geom_point(alpha = 1/3,
                 size = 0.75,
                 position = position_jitterdodge(dodge.width = 0.9)) +
      labs(x = "Gene Set Cluster",
           y = "Enrichment",
           fill = "Source") +
      theme_bw() +
      theme(
        legend.position = "bottom"
      ) +
      coord_flip()
  }
  if ( nCond == '1' ) {
    geneSetsList@cluster <- Seurat::AddMetaData(
      object = geneSetsList@cluster,
      metadata = geneSetsList@enrichmentScore$enrichment_score,
      col.name = 'enrichment_score'
    )

    df <-
      geneSetsList@cluster %>%
      #tidyseurat::mutate(label = geneSetsList@cluster@active.ident) %>%
      mutate(label = geneSetsList@cluster@active.ident) %>%
      dplyr::group_by(label) %>%
      dplyr::mutate(median_enrichment_score = median(enrichment_score)) %>%
      dplyr::ungroup() %>%
      mutate(
        label = reorder(label, enrichment_score, FUN = median)
      )

    if(removeClusterId) {
      df$label <- stringr::str_remove(df$label, "^[0-9]+: ")
    }

    p1 <-
      df %>%
      ggplot(aes(
        x = label,
        y = enrichment_score,
        fill = median_enrichment_score
      )) +
      geom_violin(scale = 'width', draw_quantiles = 0.5) +
      geom_point(alpha = 1/3,
                 size = 0.75,
                 position = position_jitterdodge(dodge.width = 0.9)) +
      labs(x = "Gene Set Cluster",
           y = "Enrichment",
           fill = "Source") +
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
            src = df$label,
            clusterNameToHighligth
          )
        ))
    )

  }


  return(p1)
}
