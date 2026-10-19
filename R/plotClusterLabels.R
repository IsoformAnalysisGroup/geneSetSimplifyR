#' Plot label information for a particular cluster
#'
#' @description Shows why a cluster got its label: the top label candidates
#'   with their combined and TF-IDF scores, the most central gene sets, the
#'   most frequent words and genes, and the largest gene sets of the cluster.
#'   Labels are the original ones from `labelClusters()`, also after a rename
#'   with `updateClusterLabels()`.
#'
#' @param geneSetsList 'gsList' object that has been clustered and labeled.
#' @param clusterNumber Which cluster to plot, given by its cluster number.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param nLabels The number of labels to plot. Default to 5
#' @param maxChar The max number of characters of labels shown in plot. Longer labels will be truncated using "...". Defaults to 50
#'
#' @return a ggplot2 object.
#' @export
#'
#' @examples
#' plotClusterLabels(geneSetsList = exampleGsList, clusterNumber = 0)

plotClusterLabels <- function(
    geneSetsList,
    clusterNumber,
    resolution = NULL,
    nLabels = 5,
    maxChar = 50
) {
  checkGsList(geneSetsList)

  resolutionColumn <- resolveResolutionColumn(geneSetsList, resolution)

  clusterIds <- levels(geneSetsList@cluster@meta.data[[resolutionColumn]])
  if (missing(clusterNumber) || length(clusterNumber) != 1L ||
      !as.character(clusterNumber) %in% clusterIds) {
    stop(
      "'clusterNumber' should be one of the clusters at resolution ",
      sub("^GSEA_snn_res\\.", "", resolutionColumn), ": ",
      paste(clusterIds, collapse = ", "), "."
    )
  }

  ### Extract plot data
  {
    ### Overwrite with original labels if necessary
    if( 'org_label' %in% colnames(geneSetsList@labels[[resolutionColumn]]$labelsDefault) ) {
      geneSetsList@labels[[resolutionColumn]]$labelsDefault$label <-
        geneSetsList@labels[[resolutionColumn]]$labelsDefault$org_label
    }

    localIntegrated <-
      geneSetsList@labels[[resolutionColumn]]$labelsDefault %>%
      dplyr::filter(cluster == clusterNumber) %>%
      dplyr::rename(value = weighted_score) %>%
      dplyr::mutate(
        label = factor(label, levels=unique(label)),
        measure = 'Label Weighted TF IDF',
        index = seq_len(dplyr::n())
      ) %>%
      dplyr::filter(index <= nLabels) %>%
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
        label = factor(label, levels=unique(label)),
        measure = 'Label TF IDF',
        index = seq_len(dplyr::n())
      ) %>%
      dplyr::filter(index <= nLabels) %>%
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
        label = factor(label, levels=unique(label)),
        measure = 'Gene Set Centrality Score',
        index = seq_len(dplyr::n())
      ) %>%
      dplyr::filter(index <= nLabels) %>%
      dplyr::select(
        label,
        value,
        measure
      )

    clusterMeta <- geneSetsList@cluster@meta.data
    gsInCluster <- clusterMeta$gene_set[
      as.character(clusterMeta[[resolutionColumn]]) == as.character(clusterNumber)
    ]

    # drop = FALSE keeps this a data.frame for a single-gene-set cluster, which
    # would otherwise collapse to a plain vector and break colSums()/apply().
    clusterMembership <-
      geneSetsList@genesetsDataframe[, gsInCluster, drop = FALSE]

    localGsSize <-
      tibble::tibble(
        label = gsInCluster,
        value = colSums(clusterMembership),
        measure = 'Gene Set Size'
      ) %>%
      dplyr::slice_max(order_by = value, n=nLabels) %>%
      dplyr::mutate(
        label = factor(label, levels=unique(label) ),
      )

    localGeneFreq <-
      apply(
        X = clusterMembership,
        MARGIN = 2,
        function(x) {
          rownames(geneSetsList@genesetsDataframe)[which(
            x == 1
          )]
        }
      ) %>%
      unlist() %>%
      table() %>%
      sort(decreasing = TRUE) %>%
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
    dplyr::filter(value != 0) %>%
    ggplot(aes(x=tidytext::reorder_within(label, value, measure), y=value)) +
    geom_bar(stat="identity", position='dodge') +
    facet_wrap(~measure, ncol=3, scales = 'free') +
    tidytext::scale_x_reordered() +
    coord_flip() +
    theme_bw() +
    labs(
      x = NULL,
      y = NULL
    )

  return(p1)
}

