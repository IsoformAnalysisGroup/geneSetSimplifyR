#' Get labels for cluster resolution
#'
#' @param geneSetsList 'gsList' object with clustered gene sets
#' @param stopWords Vector with stop words. Default is provided in 'data(stopWords)'.
#' @param clusterResolution Resolution to label.
#' @param nLabels Number of top labels to keep. Default to 10
#' @param removeFirstWord Logical. If TRUE, strip the first word (everything before
#'   the first underscore) from each gene-set name before generating bigram tokens.
#'   Useful for MSigDB-style names where the first token is an author or database
#'   prefix (e.g. "GOBP", "GOMF", "HALLMARK"). Default: FALSE.
#'
#' @return 'gsList' object with labels for selected resolution
#' @keywords internal
#' @noRd
#'
getLabels <- function(geneSetsList,
                      stopWords,
                      clusterResolution,
                      nLabels = 10,
                      removeFirstWord = FALSE
                      ) {
  seuratObj <- geneSetsList@cluster

  wordsDataset <- tibble::tibble()
  centralGenesets <- tibble::tibble()
  matrix <- as.matrix(seuratObj@graphs[["GSEA_snn"]])

  # Cluster ids are consecutive and start at 0 (see numericClusterFactor())
  clusterIds <- as.numeric(levels(seuratObj@meta.data[[clusterResolution]]))

  for (i in clusterIds + 1) {
    genesetsCluster <-
      tibble::as_tibble_col(Seurat::WhichCells(object = seuratObj, ident = i - 1),
                            column_name = "geneset")

    wordsCluster <- genesetsCluster %>%
      # Keep 'geneset' as the original name (with its underscores) so it still
      # matches the centrality table; tokenize a space-delimited copy instead.
      dplyr::mutate(tokenSource = stringr::str_replace_all(geneset, "_", " ")) %>%
      dplyr::mutate(
        tokenSource = if (removeFirstWord) {
          stringr::str_remove(tokenSource, "^\\S+\\s*")
        } else {
          tokenSource
        }
      ) %>%
      dplyr::mutate(tokenSource = stringr::str_squish(tokenSource)) %>%
      dplyr::filter(tokenSource != "") %>% # drop names left empty after stripping
      tidytext::unnest_tokens(label,
                              tokenSource,
                              token = "ngrams",
                              n = 2,
                              drop = FALSE) %>%
      # single-word names can't form a bigram (label is NA); fall back to the
      # word itself so informative one-word gene sets are still labeled.
      dplyr::mutate(label = dplyr::coalesce(label, stringr::str_to_lower(tokenSource))) %>%
      dplyr::mutate(label = stringr::str_to_title(label)) %>%
      dplyr::mutate(cluster = i - 1) %>%
      dplyr::select(-tokenSource)

    wordsDataset <- dplyr::bind_rows(wordsDataset, wordsCluster)

    clusterGenesets <-
      Seurat::WhichCells(object = seuratObj, ident = i - 1)
    # drop = FALSE: a singleton cluster would otherwise lose its gene-set name
    clusterMatrix <- matrix[clusterGenesets, clusterGenesets, drop = FALSE]

    g <-
      igraph::graph_from_adjacency_matrix(
        clusterMatrix,
        diag = FALSE,
        weighted = TRUE,
        mode = "undirected"
      )

    igraph::E(g)$distance <- 1 / igraph::E(g)$weight
    c <- igraph::betweenness(g, weights=igraph::E(g)$distance)

    # Keep the gene-set names even if the graph has no vertex names
    if (is.null(names(c))) {
      names(c) <- clusterGenesets
    }

    c <- sort(c, decreasing = TRUE)

    centralGenesetsCluster <-
      tibble::tibble(
        cluster = i - 1,
        geneset = names(c),
        centrality_score = c
      )

    centralGenesets <-
      dplyr::bind_rows(centralGenesets, centralGenesetsCluster)
  }

  geneSetsList@labels[[clusterResolution]][["labelsCentrality"]] <-
    centralGenesets %>%
    dplyr::mutate(cluster = as.numeric(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(centrality_score, n = nLabels, with_ties = FALSE) %>%
    dplyr::ungroup()

  # Drop bigrams containing a stop word; an empty list keeps every bigram
  stopWords <- stopWords[!is.na(stopWords) & nzchar(stopWords)]
  if (length(stopWords)) {
    # Lookarounds instead of \b, which fails next to non-word characters
    wordListRegex <- paste0(
      "(?<![[:alnum:]])(?:",
      paste(stringr::str_escape(stopWords), collapse = "|"),
      ")(?![[:alnum:]])"
    )
    hasStopWord <- grepl(
      wordListRegex, wordsDataset$label, ignore.case = TRUE, perl = TRUE
    )
  } else {
    hasStopWord <- rep(FALSE, nrow(wordsDataset))
  }

  cleanedWords <- wordsDataset[!hasStopWord, , drop = FALSE] %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::count(cluster, label, sort = TRUE)

  totalWords <- cleanedWords %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(total = sum(n))

  cleanedWords <-
    dplyr::left_join(cleanedWords, totalWords, by = "cluster")

  clustersTFIDF <- cleanedWords %>%
    tidytext::bind_tf_idf(label, cluster, n) %>%
    dplyr::select(cluster, label, tf_idf)

  geneSetsList@labels[[clusterResolution]][["labelsTFIDF"]] <-
    clustersTFIDF %>%
    dplyr::mutate(cluster = as.numeric(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(tf_idf, n = nLabels, with_ties = FALSE) %>%
    dplyr::ungroup()

  clustersTFIDFGenesets <- wordsDataset %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::right_join(clustersTFIDF, by = c("label", "cluster"))

  centralGenesets <- centralGenesets %>%
    dplyr::mutate(centrality_score_transformed = log10(centrality_score + 1)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      centrality_score_scaled = {
        m <- max(centrality_score_transformed, na.rm = TRUE)
        if (m > 0) centrality_score_transformed / m else 0
      }
    )%>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::ungroup()

  # Scale TF-IDF per cluster to [0, 1], the range of the centrality score
  clustersTFIDFGenesets <- clustersTFIDFGenesets %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(tf_idf_scaled = if (max(tf_idf, na.rm = TRUE) > 0) tf_idf / max(tf_idf, na.rm = TRUE) else 0) %>%
    dplyr::ungroup()

  combinedScores <-
    dplyr::inner_join(clustersTFIDFGenesets,
                      centralGenesets,
                      by = c("cluster", "geneset")) %>%
    dplyr::group_by(cluster, label, tf_idf) %>% # remove duplicate labels within clusters
    dplyr::slice(which.max(centrality_score_scaled)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(weighted_score = tf_idf_scaled + centrality_score_scaled)

  labelsDefault <- combinedScores %>%
    dplyr::mutate(cluster = as.numeric(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(weighted_score, n = nLabels, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(label, cluster, weighted_score)

  # Clusters sharing a top label get their second label appended, then their
  # id if still not unique, so every cluster name is distinct
  labelsDefault <- labelsDefault %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(weighted_score), .by_group = TRUE) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()

  topRanked <- labelsDefault %>% dplyr::filter(rank == 1)
  sharedTop <- topRanked$label[
    duplicated(topRanked$label) | duplicated(topRanked$label, fromLast = TRUE)
  ]

  if (length(sharedTop)) {
    secondLabel <- labelsDefault %>%
      dplyr::filter(rank == 2) %>%
      dplyr::select(cluster, second_label = label)

    labelsDefault <- labelsDefault %>%
      dplyr::left_join(secondLabel, by = "cluster") %>%
      dplyr::mutate(
        label = ifelse(
          rank == 1 & label %in% sharedTop & !is.na(second_label),
          paste0(label, "; ", second_label),
          label
        )
      ) %>%
      dplyr::select(-second_label)
  }

  topRanked <- labelsDefault %>% dplyr::filter(rank == 1)
  stillShared <- topRanked$label[
    duplicated(topRanked$label) | duplicated(topRanked$label, fromLast = TRUE)
  ]
  if (length(stillShared)) {
    labelsDefault <- labelsDefault %>%
      dplyr::mutate(
        label = ifelse(
          rank == 1 & label %in% stillShared,
          paste0(label, " (", cluster, ")"),
          label
        )
      )
  }

  geneSetsList@labels[[clusterResolution]][["labelsDefault"]] <-
    labelsDefault %>%
    dplyr::select(label, cluster, weighted_score)

  return(geneSetsList)
}
