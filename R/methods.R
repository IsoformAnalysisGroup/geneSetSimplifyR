#' @describeIn gsList-class Print a short summary of the object.
#' @param object A 'gsList' object.
setMethod("show", "gsList", function(object) {
  nGeneSets <- dplyr::n_distinct(object@enrichmentScore$pathway)
  nSources <- dplyr::n_distinct(object@enrichmentScore$source)
  isClustered <- length(object@cluster) > 0

  summaryLine <- paste0(
    'A gsList object with ', nGeneSets, ' gene sets from ',
    nSources, ' source(s)'
  )
  if (isClustered) {
    nClusters <- length(levels(object@cluster@active.ident))
    summaryLine <- paste0(
      summaryLine, ' clustered into ', nClusters,
      if (length(object@labels)) ' labeled clusters.' else ' clusters.'
    )
  } else {
    summaryLine <- paste0(summaryLine, '.')
  }
  cat(summaryLine, "\n", sep = "")

  meta <- if (.hasSlot(object, "metadata")) object@metadata else list()
  if (length(meta)) {
    cat("Metadata:\n")
    if (!is.null(meta$cluster_algorithm)) {
      cat("  algorithm: ", meta$cluster_algorithm, "\n", sep = "")
    }
    if (!is.null(meta$k)) {
      cat("  k (KNN): ", meta$k, "\n", sep = "")
    }
    if (!is.null(meta$lsi_elbow_cutoff)) {
      cat("  LSI elbow cutoff: ", meta$lsi_elbow_cutoff, "\n", sep = "")
    }
    if (!is.null(meta$selected_resolution)) {
      cat("  selected resolution: ", meta$selected_resolution, "\n", sep = "")
    }
    if (isClustered && !is.null(meta$active_resolution)) {
      activeResolution <- sub("^GSEA_snn_res\\.", "", meta$active_resolution)
      if (!identical(as.numeric(activeResolution),
                     as.numeric(meta$selected_resolution))) {
        cat("  active resolution: ", activeResolution, "\n", sep = "")
      }
    }
    if (!is.null(meta$n_pathways)) {
      cat("  n pathways: ", meta$n_pathways, "\n", sep = "")
    }
  }

  invisible(object)
})
