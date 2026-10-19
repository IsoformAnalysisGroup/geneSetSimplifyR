#' Coerce a cluster-id vector/factor to consecutive ids starting at 0
#'
#' Internal helper. Backends differ in how they number clusters: Seurat is
#' 0-based, igraph is 1-based, a custom backend may return gaps (e.g.
#' `c(1, 3, 5)`) or non-numeric labels. Seurat may also store factor levels in
#' lexicographic order ("10" before "2"), which breaks code that pairs labels to
#' levels by position. This maps the ids to consecutive integers starting at 0,
#' in numeric order when every id is numeric (character order otherwise), and
#' returns a factor with numeric level order.
#'
#' @param x A factor or vector of cluster ids.
#' @return A factor with levels `0, 1, ..., k - 1`.
#' @keywords internal
numericClusterFactor <- function(x) {
  x <- as.character(x)
  asNumber <- suppressWarnings(as.numeric(x))
  if (!anyNA(asNumber[!is.na(x)])) {
    uniqueIds <- sort(unique(asNumber))
    ids <- match(asNumber, uniqueIds) - 1L
  } else {
    uniqueIds <- sort(unique(x))
    ids <- match(x, uniqueIds) - 1L
  }
  factor(ids, levels = seq_along(uniqueIds) - 1L)
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
#' gsListObject <- normalizeClusterLevels(exampleGsList)
#' levels(gsListObject@cluster$GSEA_snn_res.1.1)
normalizeClusterLevels <- function(geneSetsList) {
  checkGsList(geneSetsList, needLabels = FALSE)
  resCols <- availableResolutionColumns(geneSetsList)
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
#' Three things this guards against:
#' \itemize{
#'   \item **Order.** The label table is in numeric cluster order; a factor's
#'     levels need not be, so labels are matched by cluster id rather than by
#'     position.
#'   \item **Coverage.** A cluster can be missing from the label table (e.g. a
#'     cluster whose gene-set names yield no usable bigram). Those fall back to
#'     the bare cluster id, so the returned vector always has the same length as
#'     `clusterIds`.
#'   \item **Uniqueness.** When `removeClusterId = TRUE`, two clusters with the
#'     same label (e.g. after a manual rename) would otherwise be merged into a
#'     single factor level. Those clusters keep their id prefix instead.
#' }
#'
#' @param geneSetsList A labeled 'gsList' object.
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
    dplyr::slice_max(weighted_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  clusterIds <- as.character(clusterIds)
  idx <- match(clusterIds, as.character(topLabels$cluster))
  withId <- ifelse(
    is.na(idx),
    clusterIds,
    paste0(clusterIds, ": ", topLabels$label[idx])
  )

  if (!removeClusterId) {
    return(withId)
  }

  out <- ifelse(is.na(idx), clusterIds, topLabels$label[idx])
  isDuplicated <- duplicated(out) | duplicated(out, fromLast = TRUE)
  out[isDuplicated] <- withId[isDuplicated]
  return(out)
}

#' Relabel the levels of a factor
#'
#' Internal helper replacing `plyr::mapvalues()` for factors: every level found
#' in `from` is renamed to the matching element of `to`; other levels are kept.
#'
#' @param x A factor.
#' @param from,to Character vectors of equal length.
#' @return The relabeled factor.
#' @keywords internal
relabelFactor <- function(x, from, to) {
  if (!is.factor(x)) {
    x <- factor(x)
  }
  lev <- levels(x)
  idx <- match(lev, from)
  lev[!is.na(idx)] <- to[idx[!is.na(idx)]]
  levels(x) <- lev
  return(x)
}

#' Metadata column name of a resolution
#'
#' @param resolution A single resolution, numeric or coercible to numeric.
#' @return The name of the matching column, e.g. `"GSEA_snn_res.1.1"`.
#' @keywords internal
resolutionColumnName <- function(resolution) {
  paste0("GSEA_snn_res.", as.character(as.numeric(resolution)))
}

#' Resolution columns available in a clustered gsList
#'
#' @param geneSetsList A clustered 'gsList' object.
#' @return Character vector of `GSEA_snn_res.<res>` column names.
#' @keywords internal
availableResolutionColumns <- function(geneSetsList) {
  grep(
    "^GSEA_snn_res\\.",
    names(geneSetsList@cluster@meta.data),
    value = TRUE
  )
}

#' Resolve which resolution a function should use
#'
#' Internal helper shared by the plotting and table functions. An explicit
#' `resolution` always wins. With `resolution = NULL` the active resolution is
#' used: the one last set by `updateClusterLabels()` or `updateResolution()`,
#' or the optimal resolution found by `clusterGeneSets()` if neither was used.
#'
#' @param geneSetsList A clustered 'gsList' object.
#' @param resolution `NULL` or a single resolution that was clustered.
#' @return The name of the resolution column, e.g. `"GSEA_snn_res.1.1"`.
#' @keywords internal
resolveResolutionColumn <- function(geneSetsList, resolution = NULL) {
  available <- availableResolutionColumns(geneSetsList)

  if (!is.null(resolution)) {
    resolutionNumber <- suppressWarnings(as.numeric(resolution))
    if (length(resolution) != 1L || is.na(resolutionNumber)) {
      stop("'resolution' should be a single number (or NULL).")
    }
    resolutionColumn <- resolutionColumnName(resolutionNumber)
    if (!resolutionColumn %in% available) {
      stop(
        "Resolution ", resolution, " was not used in 'clusterGeneSets()'. ",
        "Available resolution(s): ",
        paste(sub("^GSEA_snn_res\\.", "", available), collapse = ", "), "."
      )
    }
    return(resolutionColumn)
  }

  activeColumn <- NULL
  if (.hasSlot(geneSetsList, "metadata")) {
    activeColumn <- geneSetsList@metadata$active_resolution
  }
  if (!is.null(activeColumn) && activeColumn %in% available) {
    return(activeColumn)
  }

  return(geneSetsList@cluster@meta.data$best_cluster_res[1])
}

#' Check that a gsList is clustered (and optionally labeled)
#'
#' @param geneSetsList Object to check.
#' @param needLabels Whether labels from `labelClusters()` are required.
#' @return `invisible(TRUE)`; errors otherwise.
#' @keywords internal
checkGsList <- function(geneSetsList, needLabels = TRUE) {
  if (!inherits(geneSetsList, "gsList")) {
    stop("'geneSetsList' should be a gsList object.")
  }

  if (length(geneSetsList@cluster) == 0) {
    stop(
      "No clustering performed for 'geneSetsList'. ",
      "Please run 'clusterGeneSets()' and 'labelClusters()' first."
    )
  }

  if (needLabels && length(geneSetsList@labels) == 0) {
    stop(
      "No labels created for 'geneSetsList'. ",
      "Please run 'labelClusters()' first."
    )
  }

  invisible(TRUE)
}

#' Set the cluster identities used for plotting
#'
#' Internal helper. Sets `active.ident` of the Seurat object to the clusters of
#' the requested resolution, labeled as `"<id>: <label>"` (or `"<label>"` with
#' `removeClusterId = TRUE`). Only the returned copy is changed, so plotting a
#' different resolution never alters the user's object.
#'
#' @param geneSetsList A labeled 'gsList' object.
#' @param resolution `NULL` (the active resolution) or a clustered resolution.
#' @param removeClusterId Logical. Drop the `"<id>: "` prefix from the labels.
#' @return The 'gsList' object with updated identities.
#' @keywords internal
changeUsedResolution <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE
) {
  resolutionColumn <- resolveResolutionColumn(geneSetsList, resolution)

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

  geneSetsList@cluster@active.ident <- relabelFactor(
    x = geneSetsList@cluster@active.ident,
    from = currentClusterLabels,
    to = defaultLabels
  )
  return(geneSetsList)
}

#' Cluster names to highlight in a plot
#'
#' Internal helper. Translates cluster numbers into the cluster names shown in
#' the plot of the given resolution.
#'
#' @param geneSetsList A labeled 'gsList' object.
#' @param resolution `NULL` (the active resolution) or a clustered resolution.
#' @param highlightClusterNo Numeric vector of cluster numbers, or `NULL`.
#' @param removeClusterId Logical. Whether the plot drops the `"<id>: "` prefix.
#' @return Character vector of cluster names (empty when nothing to highlight).
#' @keywords internal
highlightedClusterNames <- function(
    geneSetsList,
    resolution,
    highlightClusterNo,
    removeClusterId
) {
  if (is.null(highlightClusterNo)) {
    return(character(0))
  }

  resolutionColumn <- resolveResolutionColumn(geneSetsList, resolution)
  clusterIds <- levels(geneSetsList@cluster@meta.data[[resolutionColumn]])
  highlightClusterNo <- as.character(highlightClusterNo)

  unknown <- setdiff(highlightClusterNo, clusterIds)
  if (length(unknown)) {
    warning(
      "Cluster(s) ", paste(unknown, collapse = ", "),
      " do not exist at this resolution and cannot be highlighted.",
      call. = FALSE
    )
  }

  allNames <- alignedClusterLabels(
    geneSetsList = geneSetsList,
    resolutionColumn = resolutionColumn,
    clusterIds = clusterIds,
    removeClusterId = removeClusterId
  )
  return(allNames[clusterIds %in% highlightClusterNo])
}

#' Update default resolution of gsList
#'
#' @description Sets the resolution used by all plotting and table functions
#'   when their `resolution` argument is `NULL`. It replaces the optimal
#'   resolution found by `clusterGeneSets()` as the default.
#'
#' @param geneSetsList 'gsList' object that has been clustered and labeled.
#' @param resolution Numeric. The clustering resolution to use as the new
#'   default. It must have been used as a resolution in `clusterGeneSets()`.
#'
#' @return 'gsList' object updated with the new default resolution.
#' @export
#'
#' @examples
#' gsListObject <- updateResolution(exampleGsList, 0.6)
#' extractResolution(gsListObject)
updateResolution <- function(
    geneSetsList,
    resolution
) {
  checkGsList(geneSetsList)
  geneSetsList <- upgradeGsList(geneSetsList)
  if (missing(resolution) || is.null(resolution)) {
    stop("'resolution' should be a single number.")
  }

  resolutionColumn <- resolveResolutionColumn(geneSetsList, resolution)

  geneSetsList@cluster$best_cluster_res <- resolutionColumn
  geneSetsList@metadata$selected_resolution <- as.numeric(resolution)
  geneSetsList@metadata$active_resolution <- resolutionColumn

  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE
  )

  return(geneSetsList)
}

#' Extract the resolution currently used by default
#'
#' @description Returns the resolution that the plotting and table functions
#'   use when their `resolution` argument is `NULL`. This is the optimal
#'   resolution found by `clusterGeneSets()`, unless it was changed with
#'   `updateResolution()` or `updateClusterLabels(resolution = ...)`.
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#'
#' @return The default resolution, as a number.
#' @export
#'
#' @examples
#' extractResolution(exampleGsList)
extractResolution <- function(
    geneSetsList
) {
  checkGsList(geneSetsList, needLabels = FALSE)

  res <- resolveResolutionColumn(geneSetsList, NULL) %>%
    sub(pattern = "^GSEA_snn_res\\.", replacement = "") %>%
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
#' @param geneSetsList A labeled 'gsList' object.
#' @param resolution `NULL` (the active resolution) or a clustered resolution.
#' @param scoreColumn Name of the column in the `enrichmentScore` slot to
#'   summarize, e.g. `'enrichment_score'` or `'nlog10pval'`. Should already have
#'   been checked with `validatePlotFrom()` in the caller.
#' @param summaryStat One of 'mean', 'median' or 'max'. Must already be resolved
#'   to a single value (i.e. via `match.arg()` in the caller).
#' @param removeClusterId Logical. Drop the `"<id>: "` prefix from the labels.
#' @return A data.frame with columns `source`, `cluster`, `enrichmentSummary`
#'   and `genesets` (the number of gene sets in the cluster and source,
#'   including gene sets without a value), arranged by cluster.
#' @keywords internal
summarizeEnrichment <- function(
    geneSetsList,
    resolution = NULL,
    scoreColumn = 'enrichment_score',
    summaryStat = 'mean',
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
    resolution = resolution,
    summarize = FALSE,
    removeClusterId = removeClusterId
  )

  if (! scoreColumn %in% colnames(enrichmentTable)) {
    stop("Column '", scoreColumn, "' was not found in the 'enrichmentScore' slot.")
  }

  # Copied to a fixed name so the summary below works for any requested column
  enrichmentTable$scoreToSummarize <- enrichmentTable[[scoreColumn]]

  # extractTable() returns the rows in cluster-id order; keep that order
  clusterOrder <- unique(enrichmentTable$cluster)

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
    dplyr::arrange(match(cluster, clusterOrder))
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
  ### Names ending in a direction suffix keep it after truncation
  directionPattern <- '_(up|dn|down)$'
  hasDirection <- grepl(directionPattern, aVec, ignore.case = TRUE)
  suffix <- rep('', length(aVec))
  suffix[hasDirection] <- regmatches(
    aVec[hasDirection],
    regexpr(directionPattern, aVec[hasDirection], ignore.case = TRUE)
  )

  ### Only names longer than 'maxChar' are shortened
  toModify <- nchar(aVec) > maxChar
  if (any(toModify)) {
    keepChar <- pmax(maxChar - 3 - nchar(suffix[toModify]), 1) # 3 for "..."
    aVec[toModify] <- paste0(
      substr(aVec[toModify], 1, keepChar),
      '...',
      suffix[toModify]
    )
  }

  return(aVec)
}

boldXaxisLabels <- function(src, levelsToBold) {
  if (!is.factor(src)) src <- factor(src)
  ifelse(levels(src) %in% levelsToBold, "bold", "plain")
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
#' \donttest{
#' if (requireNamespace("msigdbr", quietly = TRUE)) {
#'   # Hallmark gene sets (gene symbols)
#'   hallmark <- getGeneSets(collection = "H", gene_id_type = "gene_symbol")
#'
#'   # Just the genes of two specific pathways
#'   twoSets <- getGeneSets(
#'     pathway_names = c("HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA"),
#'     collection = "H",
#'     gene_id_type = "gene_symbol"
#'   )
#'   lengths(twoSets)
#' }
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
