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

  # Initialize variables
  wordsDataset <- tibble::tibble()
  centralGenesets <- tibble::tibble()
  matrix <- as.matrix(seuratObj@graphs[["GSEA_snn"]])

  # Cluster ids are consecutive and start at 0 (see numericClusterFactor())
  clusterIds <- as.numeric(levels(seuratObj@meta.data[[clusterResolution]]))

  # Get info for each cluster
  for (i in clusterIds + 1) {
    # GET TOKENS
    # example:
      # geneset
      # NUTT_GBM_VS_AO_GLIOMA_UP
      # KAAB_FAILED_HEART_ATRIUM_DN
    genesetsCluster <-
      tibble::as_tibble_col(Seurat::WhichCells(object = seuratObj, ident = i - 1),
                            column_name = "geneset")

    # example:
      #  cluster, label
      #  1, NUTT GBM
      #  1, GBM VS
      #  1, VS A0
      #  1, A0 GLIOMA
      # ...
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

    # GET GENESETS
    clusterGenesets <-
      Seurat::WhichCells(object = seuratObj, ident = i - 1)
    # drop = FALSE is essential: a singleton cluster would otherwise collapse to
    # a bare scalar, losing the dimnames, and the resulting one-vertex graph has
    # no vertex name -- so betweenness() returns an unnamed value, the gene set
    # is recorded as NA, and the cluster silently disappears in the join below.
    clusterMatrix <- matrix[clusterGenesets, clusterGenesets, drop = FALSE]

    # central gene set
    g <-
      igraph::graph_from_adjacency_matrix(
        clusterMatrix,
        diag = FALSE,
        weighted = TRUE,
        mode = "undirected"
      )

    igraph::E(g)$distance <- 1 / igraph::E(g)$weight
    c <- igraph::betweenness(g, weights=igraph::E(g)$distance)

    # Belt and braces: keep the gene set names attached even if the graph came
    # back without vertex names.
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

  # GET TFIDF
  # Drop bigrams containing a stop word. Words are escaped so stop words with
  # regex symbols (e.g. "c++") are matched literally; an empty stop-word list
  # keeps every bigram.
  stopWords <- stopWords[!is.na(stopWords) & nzchar(stopWords)]
  if (length(stopWords)) {
    wordListRegex <- paste0(
      "\\b(?:", paste(stringr::str_escape(stopWords), collapse = "|"), ")\\b"
    )
    hasStopWord <- grepl(wordListRegex, wordsDataset$label, ignore.case = TRUE)
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

  # SCALE CENTRALITY
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

  # SCALE TF-IDF (per cluster, to [0, 1] — matches centrality range)
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

  # DISAMBIGUATE SHARED TOP LABELS
  # The top-ranked bigram per cluster is used as the cluster's displayed name.
  # When several clusters share the same top bigram, append each cluster's
  # second-ranked bigram so the names become unique, e.g. two clusters both
  # labeled "Metabolic Process" become "Metabolic Process; Mrna Splicing" and
  # "Metabolic Process; Catabolic Process". If names are still not unique (e.g.
  # the second labels are shared too, or a cluster has no second label), the
  # cluster id is appended as a last resort, so every cluster keeps a distinct
  # name even when plots drop the "<id>: " prefix.
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
