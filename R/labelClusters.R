#' Label clusters in clustered gene sets list object
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#' @param stopWords Vector with stop words. The default vector can be viewed with data(stopWords).
#' To add new stop words 'stopWords <- c(geneSetSimplifyR::stopWords, 'new word 1', 'new word 2', ...)
#' @param nLabels Number of top labels to keep. Default to 10
#' @param removeFirstWord Logical. If TRUE, strip the first word (everything before
#'   the first underscore) from each gene-set name before generating bigram tokens.
#'   Useful for MSigDB-style names where the first token is an author or database
#'   prefix (e.g. "GOBP", "GOMF", "HALLMARK"). Default: FALSE.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @return 'gsList' object with labeled clusters for all clustered resolutions.
#' @export
#'
#' @examples
#' # Relabel the shipped example, e.g. with extra stop words
#' gsListObject <- labelClusters(
#'   geneSetsList = exampleGsList,
#'   stopWords = c(stopWords, "regulation"),
#'   verbose = FALSE
#' )
labelClusters <-
  function(geneSetsList,
           stopWords = geneSetSimplifyR::stopWords,
           nLabels = 10,
           removeFirstWord = FALSE,
           verbose = TRUE
    ) {

    # Check input
    if (!inherits(geneSetsList, "gsList")) {
      stop("'geneSetsList' should be a gsList object.")
    }
    geneSetsList <- upgradeGsList(geneSetsList)

    if (length(geneSetsList@cluster) == 0) {
      stop("No clustering performed for 'geneSetsList'. ",
           "Please run 'clusterGeneSets()' first.")
    }

    if (is.null(stopWords)) {
      stopWords <- character(0)
    }
    if (!is.character(stopWords)) {
      stop("'stopWords' should be a character vector.")
    }

    if (!is.numeric(nLabels) || length(nLabels) != 1L || nLabels < 1) {
      stop("'nLabels' should be a single positive number.")
    }

    if(verbose) {
      message("Generating labels for clusters...")
    }

    # Initialize list
    geneSetsList@labels <- list()

    resolutions <- availableResolutionColumns(geneSetsList)

    ### For each resolution extract label info
    for (i in resolutions) {
      geneSetsList@labels[[i]] <- list()

      Seurat::Idents(geneSetsList@cluster) <- geneSetsList@cluster@meta.data[[i]]

      geneSetsList <-
        getLabels(
          geneSetsList = geneSetsList,
          stopWords = stopWords,
          clusterResolution = i,
          nLabels = nLabels,
          removeFirstWord = removeFirstWord
        )

      if(verbose) {
        message("Labels created for ", i)
      }
    }

    ### Extract default labels
    defaultResolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]

    ### Reset the identities to the raw cluster ids of the default resolution.
    ### Reading them back from meta.data instead of trusting whatever
    ### active.ident happens to hold keeps labelClusters() re-runnable on an
    ### object that was already labelled.
    Seurat::Idents(geneSetsList@cluster) <-
      geneSetsList@cluster@meta.data[[defaultResolutionColumn]]
    currentClusters <- geneSetsList@cluster@active.ident
    currentClusterLabels <- levels(currentClusters)

    ### Sanity check: every labelled cluster must exist in the active identity
    labeledClusters <- as.character(unique(
      geneSetsList@labels[[defaultResolutionColumn]][["labelsDefault"]]$cluster
    ))
    if( ! all(labeledClusters %in% currentClusterLabels) ){
      stop('something went wrong with the label transfer. Contact developer')
    }

    ### Aligned by cluster id; clusters without a label keep their bare id
    ### rather than becoming NA (see alignedClusterLabels()).
    defaultLabels <- alignedClusterLabels(
      geneSetsList = geneSetsList,
      resolutionColumn = defaultResolutionColumn,
      clusterIds = currentClusterLabels,
      removeClusterId = FALSE
    )

    unlabeled <- setdiff(currentClusterLabels, labeledClusters)
    if( length(unlabeled) && verbose ) {
      message(
        "No label could be generated for cluster(s) ",
        paste(unlabeled, collapse = ", "),
        "; they keep their cluster id as label."
      )
    }

    geneSetsList@cluster@active.ident <- relabelFactor(
      x = currentClusters,
      from = currentClusterLabels,
      to = defaultLabels
    )

    ### Record labeling metadata
    geneSetsList@metadata$active_resolution <- defaultResolutionColumn
    geneSetsList@metadata$n_labels <- nLabels
    geneSetsList@metadata$labeled_at <- Sys.time()

    if(verbose) {
      message("Labeling of clusters DONE")
    }

    return(geneSetsList)
  }
