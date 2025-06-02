#' Plot label information for a particular cluster
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param clusterNumber Which cluster to plot. Here the cluster number is used
#' @param nLabels The number of labels to plot. Default to 5
#' @param maxChar The max number of characters of labels shown in plot. Longer labels will be truncated using "...". Defaults to 50
#'
#' @return a ggplot2 object.
#' @export
#'
#' @examples
#' gsListObject <- plotClusterLabels(geneSetsList = exampleGsList, clusterNumber = 0)

plotClusterLabels <- function(
    geneSetsList,
    clusterNumber,
    resolution = NULL,
    nLabels = 5,
    maxChar = 50
) {
  ### Make resolution vector
  if( ! is.null(resolution)) {
    resolutionColumn <-
      paste0("GSEA_snn_res.", resolution)
  } else {
    resolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]
  }

  ### Extract plot data
  if(TRUE) {
    ### Overwrite with original labels if nessesary
    if( 'org_label' %in% colnames(geneSetsList@labels[[resolutionColumn]]$labelsDefault) ) {
      geneSetsList@labels[[resolutionColumn]]$labelsDefault$label <-
        geneSetsList@labels[[resolutionColumn]]$labelsDefault$org_label
    }

    localIntegrated <-
      geneSetsList@labels[[resolutionColumn]]$labelsDefault %>%
      dplyr::filter(cluster == clusterNumber) %>%
      dplyr::rename(value = weighted_score) %>%
      dplyr::mutate(
        label = factor(label, levels=label),
        measure = 'Label Weighted TF IDF',
        index = 1:dplyr::n()
      ) %>%
      dplyr::filter(index %in% 1:nLabels) %>%
      dplyr::select(
        label,
        value,
        measure
      )

    localTFIDF <-
      geneSetsList@labels[[resolutionColumn]]$labelsTFIDF %>%
      dplyr::filter(cluster == clusterNumber) %>%
      dplyr::rename(value = tf_idf) %>%
      dplyr::mutate(
        label = factor(label, levels=label),
        measure = 'Label TF IDF',
        index = 1:dplyr::n()
      ) %>%
      dplyr::filter(index %in% 1:nLabels) %>%
      dplyr::select(
        label,
        value,
        measure
      )

    localCentrality <-
      geneSetsList@labels[[resolutionColumn]]$labelsCentrality %>%
      dplyr::filter(cluster == clusterNumber) %>%
      dplyr::rename(
        value = centrality_score,
        label = geneset
      ) %>%
      dplyr::mutate(
        label = factor(label, levels=label),
        measure = 'Gene Set Centrality Score',
        index = 1:dplyr::n()
      ) %>%
      dplyr::filter(index %in% 1:nLabels) %>%
      dplyr::select(
        label,
        value,
        measure
      )

    gsInCluster <-
      geneSetsList@cluster@meta.data %>%
      dplyr::filter(
        get(resolutionColumn) == clusterNumber
      ) %>%
      dplyr::pull(gene_set)

    localGsSize <-
      tibble::tibble(
        label = gsInCluster,
        value = colSums(  geneSetsList@genesetsDataframe[,gsInCluster]),
        measure = 'Gene Set Size'
      ) %>%
      dplyr::slice_max(order_by = value, n=nLabels) %>%
      dplyr::mutate(
        label = factor(label, levels=unique(label) ),
      )

    localGeneFreq <-
      apply(
        X = geneSetsList@genesetsDataframe[,gsInCluster],
        MARGIN = 2,
        function(x) {
          rownames(geneSetsList@genesetsDataframe)[which(
            x == 1
          )]
        }
      ) %>%
      unlist() %>%
      table() %>%
      sort(decreasing = T) %>%
      head(n=nLabels)
    localGeneFreqDf <-
      tibble::tibble(
        label = names(localGeneFreq),
        value = as.numeric(localGeneFreq),
        measure = 'Gene Frequency'
      ) %>%
      dplyr::mutate(
        label = factor(label, levels=unique(label) ),
      )

    topWords <-
      stringr::str_replace_all(gsInCluster, "_", " ") %>%
      stringr::str_split(pattern = ' ') %>%
      unlist() %>%
      tolower() %>%
      tibble::tibble(words = .) %>%
      dplyr::anti_join(
        x = .,
        y = tibble::tibble(words =geneSetSimplifyR::stopWords),
        by = 'words'
      ) %>%
      dplyr::group_by(words) %>%
      dplyr::summarise(
        value = dplyr::n(),
        .groups = 'drop'
      ) %>%
      dplyr::slice_max(order_by = value, n=nLabels) %>%
      dplyr::rename(label = words) %>%
      dplyr::mutate(
        measure = 'Word Frequency',
        label = factor(label, levels=unique(label) ),
      )


    jointData <-
      dplyr::bind_rows(
        localIntegrated,
        localTFIDF,
        localCentrality,
        localGsSize,
        localGeneFreqDf,
        topWords
      ) %>%
      dplyr::mutate(
        measure = factor(
          measure,
          levels = c(
            'Label Weighted TF IDF',
            'Word Frequency',
            'Gene Set Centrality Score',
            'Label TF IDF',
            'Gene Frequency',
            'Gene Set Size'
          )
        )
      )

    ### Cut to long gene names while preserving factor order
    factorOrder <- match(
      levels(jointData$label), jointData$label
    )
    jointData$label <- cutGeneSetName(
      aVec = as.character( jointData$label ),
      maxChar=maxChar
    )
    jointData$label <- factor(
      jointData$label,
      levels = unique(jointData$label[factorOrder])
    )

  }

  ### Plot
  p1 <-
    jointData %>%
    ggplot(aes(x=label, y=value)) +
    geom_bar(stat="identity", position='dodge') +
    facet_wrap(~measure, ncol=3, scales = 'free') +
    coord_flip() +
    theme_bw() +
    labs(
      x = NULL,
      y = NULL
    )

  return(p1)
}

