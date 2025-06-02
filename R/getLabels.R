#' Get labels for cluster resolution
#'
#' @param geneSetsList 'gsList' object with clustered gene sets
#' @param stopWords Vector with stop words. Default is provided in 'data(stopWords)'.
#' @param clusterResolution Resolution to label.
#' @param nLabels Number of top labels to keep. Default to 10
#'
#' @return 'gsList' object with labels for selected resolution
#' @keywords internal
#' @noRd
#'
#' @examples
#' gsListObject <- getLabels(geneSetsList = exampleGsList, stopWords = stopWords, clusterResolution = 2.1)
getLabels <- function(geneSetsList,
                      stopWords,
                      clusterResolution,
                      nLabels = 10
                      ) {
  seuratObj <- geneSetsList@cluster

  # Initialize variables
  wordsDataset <- tibble::tibble()
  centralGenesets <- tibble::tibble()
  matrix <- as.matrix(seuratObj@graphs[["GSEA_snn"]])

  # Get info for each cluster
  for (i in 1:length(levels(seuratObj@meta.data[[clusterResolution]]))) {
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
      dplyr::mutate(geneset = stringr::str_replace_all(geneset, "_", " ")) %>%
      tidytext::unnest_tokens(label,
                              geneset,
                              token = "ngrams",
                              n = 2,
                              drop = F) %>%
      dplyr::mutate(cluster = i - 1)

    wordsDataset <- dplyr::bind_rows(wordsDataset, wordsCluster)

    # GET GENESETS
    clusterGenesets <-
      Seurat::WhichCells(object = seuratObj, ident = i - 1)
    clusterMatrix <- matrix[clusterGenesets, clusterGenesets]

    # central gene set
    g <-
      igraph::graph.adjacency(
        clusterMatrix,
        diag = F,
        weighted = T,
        mode = "undirected"
      )
    c <- igraph::betweenness(g) %>%
      sort(decreasing = T)

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
    dplyr::slice_max(centrality_score, n = nLabels, with_ties = F) %>%
    dplyr::ungroup()

  # GET TFIDF
  wordListRegex = paste0("\\b(?:", paste(stopWords, collapse = "|"), ")\\b")

  cleanedWords <- wordsDataset %>%
    dplyr::filter(!grepl(wordListRegex, label)) %>%
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
    dplyr::slice_max(tf_idf, n = nLabels, with_ties = F) %>%
    dplyr::ungroup()

  clustersTFIDFGenesets <- wordsDataset %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::mutate(geneset = stringr::str_replace_all(geneset, " ", "_")) %>%
    dplyr::right_join(clustersTFIDF, by = c("label", "cluster"))

  # SCALE CENTRALITY
  centralGenesets <- centralGenesets %>%
    dplyr::mutate(centrality_score_transformed = log10(centrality_score + 1)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(centrality_score_scaled = centrality_score_transformed + 1 / max(centrality_score_transformed + 1, na.rm = TRUE)) %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::ungroup()

  combinedScores <-
    dplyr::inner_join(clustersTFIDFGenesets,
                      centralGenesets,
                      by = c("cluster", "geneset")) %>%
    dplyr::group_by(cluster, label, tf_idf) %>% # remove duplicate labels within clusters
    dplyr::slice(which.max(centrality_score_scaled)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(weighted_score = tf_idf * centrality_score_scaled)

  geneSetsList@labels[[clusterResolution]][["labelsDefault"]] <-
    combinedScores %>%
    dplyr::mutate(cluster = as.numeric(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(weighted_score, n = nLabels, with_ties = F) %>%
    dplyr::ungroup() %>%
    dplyr::select(label, cluster, weighted_score)

  return(geneSetsList)
}
