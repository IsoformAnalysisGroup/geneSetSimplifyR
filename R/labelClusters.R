#' Label clusters in clustered gene sets list object
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#' @param stopWords Vector with stop words. The default vector can be viewed with data(stopWords).
#' To add new stop words 'stopWords <- c(geneSetSimplifyR::stopWords, 'new word 1', 'new word 2', ...)
#' @param nLabels Number of top labels to keep. Default to 10
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @return 'gsList' object with labelled clusters for all clustered resolutions.
#' @export
#'
#' @examples
#' gsListObject <- labelClusters(geneSetsList = gsListObject)
labelClusters <-
  function(geneSetsList,
           removeClusterId = FALSE,
           stopWords = geneSetSimplifyR::stopWords,
           nLabels = 10,
           verbose = TRUE
    ) {

    # Check input
    if (!class(geneSetsList) == "gsList") {
      stop("'geneSetsList' should be a gsList object.")
    }

    if (length(geneSetsList@cluster) == 0) {
      stop("No clustering performed for 'geneSetsList'.
         Please run 'clusterGeneSets()' and 'labelClusters()' before plotting.")
    }

    if (!is.logical(removeClusterId)) {
      stop("'removeClusterId' should be logical.")
    }

    if (!is.vector(stopWords)) {
      stop("'stopWords' should be a vector.")
    }

    if(verbose) {
      message("Generating labels for clusters...")
    }

    # Initialize list
    geneSetsList@labels <- list()

    resolutions <-
      grep("GSEA_snn_res.",
           names(geneSetsList@cluster@meta.data),
           value = TRUE)

    ### Extract current values
    currentClusters <-  geneSetsList@cluster@active.ident
    currentClusterLabels <-
      levels(geneSetsList@cluster@active.ident)

    ### For each resolution extract label info
    for (i in resolutions) {
      geneSetsList@labels[[i]] <- list()

      Seurat::Idents(geneSetsList@cluster) <- geneSetsList@cluster@meta.data[[i]]

      geneSetsList <-
        getLabels(
          geneSetsList = geneSetsList,
          stopWords = stopWords,
          clusterResolution = i,
          nLabels = nLabels
        )

      if(verbose) {
        message(paste0("Labels created for ", i))
      }
    }

    ### Extract default labels
    defaultResolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]
    defaultLabels <-
      geneSetsList@labels[[defaultResolutionColumn]][["labelsDefault"]] %>%
      dplyr::mutate(cluster = as.numeric(cluster)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(weighted_score, n = 1, with_ties = F) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(label = paste0(cluster, ": ", label)) %>%
      dplyr::pull(label)

    ### Sanity check
    labelStrip <- sapply(
      stringr::str_split(defaultLabels, pattern = ': '),
      function(x) x[1]
    )
    if( ! all(labelStrip %in% currentClusterLabels) ){
      stop('something went wrong with the label transfer. Contact developer')
    }
    defaultLabels <- defaultLabels[match(
      currentClusterLabels, labelStrip
    )]

    geneSetsList@cluster@active.ident <-
      plyr::mapvalues(x = currentClusters,
                      from = currentClusterLabels,
                      to = defaultLabels)

    if(verbose) {
      message("Labeling of clusters DONE")
    }

    return(geneSetsList)
  }
