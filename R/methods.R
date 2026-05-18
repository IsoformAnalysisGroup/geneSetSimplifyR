setMethod("show", "gsList", function(object) {
  if( length( object@labels ) ) {
    print(
      stringr::str_c(
        'A gsList object with ',
        dplyr::n_distinct(object@enrichmentScore$pathway),
        ' gene sets from ',
        dplyr::n_distinct(object@enrichmentScore$source),
        ' source(s)',
        ' clustered into ',
        length(levels( object@cluster@active.ident )),
        ' labeled clusters.'
      )
    )
  } else if( !is.null(object@cluster) ) {
    print(
      stringr::str_c(
        'A gsList object with ',
        dplyr::n_distinct(object@enrichmentScore$pathway),
        ' gene sets from ',
        dplyr::n_distinct(object@enrichmentScore$source),
        ' source(s)',
        ' clustered into ',
        length(levels( object@cluster@active.ident )),
        ' clusters.'
      )
    )
  } else {
    print(
      stringr::str_c(
        'A gsList object with ',
        dplyr::n_distinct(object@enrichmentScore$pathway),
        ' gene sets from ',
        dplyr::n_distinct(object@enrichmentScore$source),
        ' source(s)'
      )
    )
  }

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
    if (!is.null(meta$n_pathways)) {
      cat("  n pathways: ", meta$n_pathways, "\n", sep = "")
    }
  }

})
