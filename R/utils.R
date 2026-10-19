#' Coerce a cluster-id vector/factor to start at 0 with numeric level order
#'
#' Internal helper. Seurat may store cluster factor levels in lexicographic
#' order ("10" before "2"), which breaks code that pairs labels to levels by
#' position. This renumbers ids to start at 0 and rebuilds the factor with
#' numeric level order.
#'
#' @param x A factor or vector of (coercible-to-numeric) cluster ids.
#' @return A factor with levels in numeric order, starting at 0.
#' @keywords internal
numericClusterFactor <- function(x) {
  ids <- as.numeric(as.character(x))
  ids <- ids - min(ids)                       # start from 0 (robust to 0/1-indexing)
  factor(ids, levels = sort(unique(ids)))     # numeric level order
}

#' Normalize cluster factor level order in a clustered gsList
#'
#' @description Rebuilds every `GSEA_snn_res.<res>` metadata column so its
#'   factor levels are in numeric order (0, 1, 2, ...), starting at 0. Newer
#'   objects get this automatically inside `clusterGeneSets()`; use this to
#'   migrate objects clustered before that fix without re-clustering. Cluster
#'   assignments are unchanged -- only the level order is corrected.
#'
#' @param geneSetsList A clustered 'gsList' object.
#' @return The 'gsList' with normalized cluster factor level order.
#' @export
#'
#' @examples
#' \dontrun{
#' neftelGss <- normalizeClusterLevels(neftelGss)
#' }
normalizeClusterLevels <- function(geneSetsList) {
  if (!inherits(geneSetsList, "gsList")) {
    stop("'geneSetsList' should be a gsList object")
  }
  resCols <- grep(
    "GSEA_snn_res.",
    names(geneSetsList@cluster@meta.data),
    value = TRUE
  )
  for (col in resCols) {
    geneSetsList@cluster@meta.data[[col]] <-
      numericClusterFactor(geneSetsList@cluster@meta.data[[col]])
  }
  return(geneSetsList)
}

#' Build display labels aligned to a set of cluster ids
#'
#' Internal helper. Takes the `labelsDefault` table for one resolution, keeps
#' the top-scoring label per cluster, and returns one label per element of
#' `clusterIds`, in that order.
#'
#' Two things this guards against:
#' \itemize{
#'   \item **Order.** The label table is in numeric cluster order; a factor's
#'     levels need not be, so labels are matched by cluster id rather than by
#'     position.
#'   \item **Coverage.** A cluster can be missing from the label table (e.g. a
#'     cluster whose gene-set names yield no usable bigram). Those fall back to
#'     the bare cluster id, so the returned vector always has the same length as
#'     `clusterIds` instead of blowing up downstream in `plyr::mapvalues()`.
#' }
#'
#' @param geneSetsList A labelled 'gsList' object.
#' @param resolutionColumn Name of the resolution, e.g. `"GSEA_snn_res.1"`.
#' @param clusterIds Character vector of cluster ids to label, in the desired
#'   output order (typically `levels(active.ident)`).
#' @param removeClusterId Logical. Drop the `"<id>: "` prefix from the labels.
#' @return Character vector of labels, one per element of `clusterIds`.
#' @keywords internal
alignedClusterLabels <- function(
    geneSetsList,
    resolutionColumn,
    clusterIds,
    removeClusterId = FALSE
) {
  labelTable <- geneSetsList@labels[[resolutionColumn]][["labelsDefault"]]

  if (is.null(labelTable)) {
    stop(
      "No labels found for '", resolutionColumn,
      "'. Please run 'labelClusters()' before plotting."
    )
  }

  topLabels <- labelTable %>%
    dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(weighted_score, n = 1, with_ties = F) %>%
    dplyr::ungroup()

  clusterIds <- as.character(clusterIds)
  idx <- match(clusterIds, as.character(topLabels$cluster))

  if (removeClusterId) {
    out <- ifelse(is.na(idx), clusterIds, topLabels$label[idx])
  } else {
    out <- ifelse(
      is.na(idx),
      clusterIds,
      paste0(clusterIds, ": ", topLabels$label[idx])
    )
  }

  return(out)
}

changeUsedResolution <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE
) {
  if( ! is.null(resolution)) {
    resolutionColumn <-
      paste0("GSEA_snn_res.", resolution)
  } else {
    resolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]
  }


  Seurat::Idents(geneSetsList@cluster) <-
    geneSetsList@cluster@meta.data[[resolutionColumn]]

  currentClusterLabels <-
    levels(geneSetsList@cluster@active.ident)

  # Match labels to the active.ident levels by cluster id (see
  # alignedClusterLabels()) -- a positional pairing breaks whenever a cluster
  # has no label or the level order is not numeric.
  defaultLabels <- alignedClusterLabels(
    geneSetsList = geneSetsList,
    resolutionColumn = resolutionColumn,
    clusterIds = currentClusterLabels,
    removeClusterId = removeClusterId
  )

  geneSetsList@cluster@active.ident <-
    plyr::mapvalues(x = geneSetsList@cluster@active.ident,
                    from = currentClusterLabels,
                    to = defaultLabels)
  return(geneSetsList)
}

#' Update default resolution of gsList
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#' @param resolution Float. The clustering resolution to use as the new default.
#' The new 'resolution' must have been used as a resolution in the 'clusterGeneSets()' function.
#'
#' @return 'gsList' object updated with new default resolution
#' @export
#'
#' @examples
#' gsListObject <- updateResolution(exampleGsList, 1.1)
updateResolution <- function(
    geneSetsList,
    resolution
) {
  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = resolution,
    removeClusterId = FALSE
  )

  geneSetsList@cluster$best_cluster_res <-
    stringr::str_c('GSEA_snn_res.',resolution)

  geneSetsList@metadata$selected_resolution <- resolution

  return(geneSetsList)
}

#' Extract default resolution
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#'
#' @return default resolution
#' @export
#'
#' @examples
#' gsListObject <- updateResolution(exampleGsList, 1.1)
extractResolution <- function(
    geneSetsList
) {

  # Check input
  if (!inherits(geneSetsList, "gsList")) {
    stop("'geneSetsList' should be a gsList object")
  }

  if (length(geneSetsList@cluster) == 0) {
    stop(
      "No clustering performed for 'geneSetsList'.
         Please run 'clusterGeneSets()' and 'labelClusters()' before plotting."
    )
  }

  res <- geneSetsList@cluster$best_cluster_res[1] %>%
    stringr::str_remove('GSEA_snn_res.') %>%
    as.numeric()

  return(res)
}


#' Validate the columns requested for plotting
#'
#' Internal helper. Used by the plotting functions to check the `plotFrom`
#' argument before anything is summarized or drawn. A column is only usable if
#' it is present in the `enrichmentScore` slot, holds at least one non-NA value,
#' and is numeric.
#'
#' Only one column is plotted. Several can be supplied - the plot functions
#' declare their candidates as a default vector, the same way `summaryStat`
#' does - in which case the first usable one is returned. That way a default of
#' `c('nlog10pval', 'enrichment_score')` still works on a `gsList` built without
#' p-values, while a single explicitly named column that cannot be plotted
#' errors rather than silently showing something else.
#'
#' @param geneSetsList A 'gsList' object.
#' @param plotFrom Name(s) of the column(s) in the `enrichmentScore` slot to
#'   plot, e.g. `'enrichment_score'`.
#' @return A single column name.
#' @keywords internal
validatePlotFrom <- function(geneSetsList, plotFrom) {
  if (!is.character(plotFrom) || !length(plotFrom) || anyNA(plotFrom)) {
    stop(
      "'plotFrom' should be a character vector with the name of a column in ",
      "the 'enrichmentScore' slot."
    )
  }

  plotFrom <- unique(plotFrom)

  enrichmentDF <- geneSetsList@enrichmentScore

  problems <- vapply(
    plotFrom,
    function(column) plotFromProblem(enrichmentDF, column),
    character(1)
  )

  ### Use the first column that can actually be plotted
  isUsable <- is.na(problems)
  if (any(isUsable)) {
    return(plotFrom[which(isUsable)[1]])
  }

  ### Nothing left - report why each candidate was rejected
  availableColumns <- setdiff(colnames(enrichmentDF), c('source', 'pathway'))

  stop(paste0(
    'Cannot plot ', quoteAndCollapse(plotFrom), ':\n',
    paste0("  '", plotFrom, "' - ", problems, collapse = '\n'), '\n',
    '  Available column(s): ', quoteAndCollapse(availableColumns), '.\n',
    "  Other scores need to be added via the 'geneSetsDF' argument of ",
    "'initializeList()' (which means redoing the workflow)."
  ))
}

#' Why a column cannot be plotted
#'
#' Internal helper for `validatePlotFrom()`.
#'
#' The three checks are done in this order on purpose: a score that was never
#' supplied to `initializeList()` is stored as a logical `NA` column, so
#' checking for "all NA" before "is numeric" gives the informative reason rather
#' than a confusing complaint about the column type.
#'
#' @param enrichmentDF The `enrichmentScore` slot of a 'gsList' object.
#' @param column Name of a single column.
#' @return `NA_character_` if the column can be plotted, otherwise a one-line
#'   description of the problem.
#' @keywords internal
plotFromProblem <- function(enrichmentDF, column) {
  if (! column %in% colnames(enrichmentDF)) {
    return("not found in the 'enrichmentScore' slot.")
  }

  if (all(is.na(enrichmentDF[[column]]))) {
    return('all values are NA.')
  }

  if (! is.numeric(enrichmentDF[[column]])) {
    return('only numeric columns can be plotted.')
  }

  return(NA_character_)
}

#' Quote and comma-collapse a character vector for use in messages
#'
#' Internal helper.
#'
#' @param x A character vector.
#' @return A single string, e.g. `"'pval', 'padj'"`. Returns `"none"` for an
#'   empty vector so error messages stay readable.
#' @keywords internal
quoteAndCollapse <- function(x) {
  if (!length(x)) {
    return('none')
  }
  paste0("'", x, "'", collapse = ", ")
}

#' Summarize enrichment per cluster and source
#'
#' Internal helper. Summarizes the per-gene-set enrichment scores of
#' `extractTable()` into one value per cluster and source, using the statistic
#' requested by the plot function. Done here rather than reusing the
#' `enrichment mean` column of `extractTable()` so the summary can follow
#' `summaryStat`.
#'
#' @param geneSetsList A labelled 'gsList' object, with the resolution to plot
#'   already applied (see `changeUsedResolution()`).
#' @param scoreColumn Name of the column in the `enrichmentScore` slot to
#'   summarize, e.g. `'enrichment_score'` or `'nlog10pval'`. Should already have
#'   been checked with `validatePlotFrom()` in the caller.
#' @param summaryStat One of 'mean', 'median' or 'max'. Must already be resolved
#'   to a single value (i.e. via `match.arg()` in the caller).
#' @param removeClusterId Logical. Drop the `"<id>: "` prefix from the labels.
#' @return A data.frame with columns `source`, `cluster`, `enrichmentSummary`
#'   and `genesets`, arranged by cluster.
#' @keywords internal
summarizeEnrichment <- function(
    geneSetsList,
    scoreColumn = 'enrichment_score',
    summaryStat = 'median',
    removeClusterId = FALSE
) {
  summaryFun <- switch(
    summaryStat,
    mean   = mean,
    median = median,
    max    = max,
    stop("Unknown 'summaryStat': ", summaryStat)
  )

  enrichmentTable <- extractTable(
    geneSetsList = geneSetsList,
    summarize = FALSE,
    removeClusterId = removeClusterId
  )

  if (! scoreColumn %in% colnames(enrichmentTable)) {
    stop("Column '", scoreColumn, "' was not found in the 'enrichmentScore' slot.")
  }

  # Copied to a fixed name so the summary below works for any requested column
  enrichmentTable$scoreToSummarize <- enrichmentTable[[scoreColumn]]

  enrichmentTable %>%
    dplyr::group_by(source, cluster) %>%
    dplyr::summarise(
      enrichmentSummary = round(
        summarizeScores(scoreToSummarize, summaryFun),
        digits = 2
      ),
      genesets = dplyr::n(),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(cluster)
}

#' Apply a summary statistic while tolerating missing values
#'
#' Internal helper. A column only has to be *partly* annotated to pass
#' `validatePlotFrom()`, so NAs are dropped before summarizing rather than
#' propagating to the whole cluster. A cluster/source with no annotated gene set
#' at all gives `NA` instead of the `-Inf`/warning `max(numeric(0))` would give.
#'
#' @param x Numeric vector.
#' @param summaryFun Function to apply, e.g. `mean`.
#' @return A single number, or `NA_real_` if `x` holds no non-NA value.
#' @keywords internal
summarizeScores <- function(x, summaryFun) {
  x <- x[!is.na(x)]

  if (!length(x)) {
    return(NA_real_)
  }

  summaryFun(x)
}

#' Axis/legend text describing an enrichment summary
#'
#' Internal helper.
#'
#' @param summaryStat One of 'mean', 'median' or 'max'.
#' @param scoreColumn Name of the summarized column. `'enrichment_score'` is
#'   spelled out as "Enrichment"; any other column is named as-is so the axis
#'   says what is actually plotted.
#' @return A single string.
#' @keywords internal
enrichmentSummaryLabel <- function(
    summaryStat = 'mean',
    scoreColumn = 'enrichment_score'
) {
  statLabel <- switch(
    summaryStat,
    mean   = 'Mean',
    median = 'Median',
    max    = 'Max',
    ''
  )

  scoreLabel <- switch(
    scoreColumn,
    enrichment_score = 'Enrichment',
    nlog10pval       = '-log10(p-value)',
    scoreColumn
  )

  trimws(paste(statLabel, scoreLabel))
}


countConditions <- function(
    geneSetsList
) {
  nCond <- 'Not assigned'

  condCount <- dplyr::n_distinct(geneSetsList@enrichmentScore$source)

  if( condCount >= 3 ) {
    nCond <- '3+'
  }
  if( condCount == 1 ) {
    nCond <- '1'
  }
  if( condCount == 2 ) {
    nCond <- '2'
  }
  return(nCond)
}

cutGeneSetName <- function(aVec, maxChar = 60) {
  # aVec <- c(paste0(rep('a', 30), collapse = ''), paste0(rep('b', 50), collapse = ''), paste0(paste0(rep('c', 50), collapse = ''), '_up'))

  ### make vector with max lengths
  maxVecCar <- rep(maxChar, length(aVec)) - 3 # for the "..."

  ### Annotat up and down
  hasDirection <- which( grepl('_up$|_dn$|_down$', aVec ))
  if( length(hasDirection) ) {
    orgDirection <- sapply(
      strsplit(aVec[hasDirection], '_' ),
      function(x) tail(x,1)
    )

    # Update max lengths
    maxVecCar[hasDirection] <- maxVecCar[hasDirection] - (nchar(orgDirection) +1) # +1 for the "_"
  }


  ### Modify those which are to long
  toModify <- which(nchar(aVec) > maxVecCar)
  if(length(toModify)) {
    aVec[toModify] <- paste0(
      substr(aVec[toModify], 1, maxVecCar[toModify]),
      '...'
    )
  }

  ### Add up, down back
  addDirBack <- intersect(hasDirection, toModify)
  if( length(addDirBack) ) {
    aVec[addDirBack] <- paste0(
      aVec[addDirBack],
      '_',
      orgDirection
    )
  }

  return(aVec)
}

boldXaxisLabels <- function(src, levelsToBold) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- levelsToBold %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(levelsToBold, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'levelsToBold' must be in src")
  }
}


#' Fetch MSigDB gene sets by pathway name
#'
#' @description Fetches gene-set memberships from MSigDB via
#'   [msigdbr::msigdbr()] and returns them in the named-list format expected
#'   by `geneSetSimplifyR()` (pathway name -> character vector of gene IDs).
#'
#'   Useful in two situations:
#'   \itemize{
#'     \item Inspecting the gene members of specific pathways (e.g. the top
#'       pathways from a cluster after running the workflow).
#'     \item Building a `geneSetsList` from scratch without going through
#'       `pairedGSEA::prepare_msigdb()`.
#'   }
#'
#' @param pathway_names Character vector of MSigDB pathway names
#'   (i.e. `gs_name` values). If `NULL` (default), all pathways in the
#'   requested collection/subcollection are returned.
#' @param species Species, passed to `msigdbr::msigdbr()`. Default
#'   `"Homo sapiens"`.
#' @param collection MSigDB collection (e.g. `"C5"`, `"H"`, `"C2"`). Default
#'   `"C5"`.
#' @param subcollection Optional MSigDB sub-collection (e.g. `"GO:BP"`).
#'   `NULL` (default) returns the full collection.
#' @param gene_id_type Which gene identifier to use for the returned vectors.
#'   One of `"ensembl_gene"` (default), `"gene_symbol"`, or `"entrez_gene"`.
#'
#' @return A named list where each element is a character vector of gene IDs
#'   for a pathway. Names are pathway names (`gs_name`).
#'
#' @examples
#' \dontrun{
#' # All pathways in C5 / GO:BP
#' allBP <- getGeneSets(collection = "C5", subcollection = "GO:BP")
#'
#' # Just the genes for two specific pathways
#' subset <- getGeneSets(
#'   pathway_names = c("GOBP_CELL_CYCLE", "GOBP_DNA_REPAIR")
#' )
#' }
#'
#' @export
getGeneSets <- function(
    pathway_names = NULL,
    species       = "Homo sapiens",
    collection    = "C5",
    subcollection = NULL,
    gene_id_type  = c("ensembl_gene", "gene_symbol", "entrez_gene")
) {
  gene_id_type <- match.arg(gene_id_type)

  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "Package 'msigdbr' is required for getGeneSets(). Install it with ",
      "install.packages(\"msigdbr\")."
    )
  }

  tib <- msigdbr::msigdbr(
    species       = species,
    collection    = collection,
    subcollection = subcollection
  )

  if (!nrow(tib)) {
    stop(sprintf(
      "No gene sets returned for species='%s', collection='%s', subcollection=%s.",
      species, collection,
      if (is.null(subcollection)) "NULL" else sprintf("'%s'", subcollection)
    ))
  }

  if (!gene_id_type %in% colnames(tib)) {
    stop(sprintf(
      "Column '%s' not found in msigdbr output. Available columns include: %s.",
      gene_id_type,
      paste(intersect(colnames(tib),
                      c("ensembl_gene", "gene_symbol", "entrez_gene")),
            collapse = ", ")
    ))
  }

  if (!is.null(pathway_names)) {
    keep <- tib$gs_name %in% pathway_names
    if (!any(keep)) {
      stop(
        "None of the requested 'pathway_names' were found in the specified ",
        "collection/subcollection."
      )
    }
    missing <- setdiff(pathway_names, tib$gs_name[keep])
    if (length(missing)) {
      warning(sprintf(
        "%d pathway name(s) not found and will be skipped (first 5: %s).",
        length(missing),
        paste(utils::head(missing, 5), collapse = ", ")
      ), call. = FALSE)
    }
    tib <- tib[keep, , drop = FALSE]
  }

  split(tib[[gene_id_type]], tib$gs_name)
}
