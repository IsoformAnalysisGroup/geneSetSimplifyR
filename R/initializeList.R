#' Initialize gene sets list
#'
#' @description Initialize a 'gsList' object from a dataframe with gene sets and enrichment scores and a gene sets database.
#'
#' @param geneSetsList A list containing the gene sets to be analyzed.
#' The name of each element in the list should be a gene set and each element should contain
#' the genes belonging to that gene set.
#' @param geneSetsDF Optional data.frame with gene sets and enrichment scores that can be useful
#' for downstream visualization (both for effect size and direction and if multiple conditions are analyzed).
#' Should have a column named 'pathway' containing the gene sets and a column
#' named 'enrichment_score' or 'NES' with the enrichment of each gene set.
#' Optionally have a column named 'pval' or 'pvalue' with the p-value of each gene set; a
#' 'nlog10pval' column with the -log10 transformed p-values is added automatically and is
#' used by default by the plotting functions that summarize scores (e.g. 'plotHeatmap()').
#' p-values must lie between 0 and 1. A p-value of 0 is transformed as the
#' smallest positive double (giving about 307.7) rather than infinity.
#' Optionally have a column named 'source' indicating the source of each gene set. This allows for
#' analysis of multi-conditional data throughout the entire workflow and is also supported by many
#' plotting functions. If a gene-set was found in multiple conditions there should be multiple rows
#' for the same 'pathway' indicating its respective source and corresponding enrichment_score.
#' Each combination of 'pathway' and 'source' may appear only once.
#' Any other numeric column is kept and can be plotted via the 'plotFrom' argument of the
#' plotting functions. List columns (e.g. 'leadingEdge' from fgsea) are dropped.
#' See the vignette for examples.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @details Gene sets without any gene are removed with a warning. Missing gene
#'   ids (`NA` or empty strings) are ignored.
#'
#' @return 'gsList' object with a gene sets dataframe and enrichment scores for each gene set.
#' @export
#'
#' @examples
#' gsListObject <- initializeList(geneSetsList = exampleGeneSets, geneSetsDF = exampleEnrichment)
#' gsListObject
initializeList <- function(
    geneSetsList,
    geneSetsDF = NULL,
    verbose = TRUE
) {

  if(verbose) {
    message("Checking input...")
  }

  ### Gene sets
  if (!is.list(geneSetsList) || is.data.frame(geneSetsList) ||
      !length(geneSetsList)) {
    stop("'geneSetsList' should be a non-empty named list of gene sets.")
  }

  geneSetNames <- names(geneSetsList)
  if (is.null(geneSetNames) || anyNA(geneSetNames) ||
      any(!nzchar(geneSetNames))) {
    stop(
      "Every element of 'geneSetsList' should be named with its gene set ",
      "name."
    )
  }

  duplicatedNames <- unique(geneSetNames[duplicated(geneSetNames)])
  if (length(duplicatedNames)) {
    stop(
      "Gene set names in 'geneSetsList' should be unique. Duplicated: ",
      quoteAndCollapse(utils::head(duplicatedNames, 5)),
      if (length(duplicatedNames) > 5) ", ..." else "", "."
    )
  }

  if (!all(vapply(geneSetsList, is.character, logical(1)))) {
    stop("each element in the 'geneSetsList' should contain a character vector")
  }

  ### Drop missing gene ids and duplicated genes within a set
  geneSetsList <- lapply(geneSetsList, function(genes) {
    unique(genes[!is.na(genes) & nzchar(genes)])
  })

  ### Enrichment table
  no_table <- FALSE
  if (is.null(geneSetsDF)) {
    geneSetsDF <- data.frame(
      pathway = names(geneSetsList),
      enrichment_score = NA,
      stringsAsFactors = FALSE
    )
    no_table <- TRUE
  }

  if (!is.data.frame(geneSetsDF)) {
    stop("'geneSetsDF' should be a data.frame (or NULL).")
  }

  # Plain data.frame: avoids data.table/tibble specific subsetting rules
  geneSetsDF <- as.data.frame(geneSetsDF, stringsAsFactors = FALSE)

  # Check pathway
  if(! 'pathway' %in% colnames(geneSetsDF)) {
    stop("Column 'pathway' was not found in geneSetsDF.")
  }
  geneSetsDF$pathway <- as.character(geneSetsDF$pathway)

  ### Drop list columns such as fgsea's 'leadingEdge'
  isListColumn <- vapply(geneSetsDF, is.list, logical(1))
  if (any(isListColumn)) {
    if (verbose) {
      message(
        "Dropping list column(s) from 'geneSetsDF': ",
        quoteAndCollapse(colnames(geneSetsDF)[isListColumn]), "."
      )
    }
    geneSetsDF <- geneSetsDF[, !isListColumn, drop = FALSE]
  }

  # make unique
  geneSetsDF <- unique(geneSetsDF)

  # Check enrichment_score
  if(! 'enrichment_score' %in% colnames(geneSetsDF)) {
    if( 'NES' %in% colnames(geneSetsDF) ) {
      if (verbose) {
        message("enrichment_score was not found in the dataframe. Use NES instead.")
      }
      geneSetsDF$enrichment_score <- geneSetsDF$NES
      geneSetsDF$NES <- NULL
    } else {
      geneSetsDF$enrichment_score <- NA
    }
  }
  if (!all(is.na(geneSetsDF$enrichment_score)) &&
      !is.numeric(geneSetsDF$enrichment_score)) {
    stop("Column 'enrichment_score' of 'geneSetsDF' should be numeric.")
  }

  # Check pvalues
  if(! 'pval' %in% colnames(geneSetsDF)) {
    if( 'pvalue' %in% colnames(geneSetsDF) ) {
      if (verbose) {
        message("pval was not found in the dataframe. Use pvalue instead.")
      }
      geneSetsDF$pval <- geneSetsDF$pvalue
      geneSetsDF$pvalue <- NULL
    } else {
      geneSetsDF$pval <- NA
    }
  }

  if (!all(is.na(geneSetsDF$pval))) {
    if (!is.numeric(geneSetsDF$pval)) {
      stop("Column 'pval' (or 'pvalue') of 'geneSetsDF' should be numeric.")
    }
    if (any(geneSetsDF$pval < 0 | geneSetsDF$pval > 1, na.rm = TRUE)) {
      stop("Column 'pval' (or 'pvalue') of 'geneSetsDF' should be between 0 and 1.")
    }
  }

  # Add the -log10 transformed p-values used by the plotting functions.
  # p-values of 0 are floored at the smallest representable double instead of
  # becoming Inf; an absent p-value column gives an all-NA column, which the
  # plotting functions skip when choosing what to plot.
  if(! 'nlog10pval' %in% colnames(geneSetsDF)) {
    geneSetsDF$nlog10pval <- -log10(
      pmax(as.numeric(geneSetsDF$pval), .Machine$double.xmin)
    )
  }

  ### Source
  if( 'source' %in% colnames(geneSetsDF) ) {
    if (anyNA(geneSetsDF$source)) {
      stop("Column 'source' of 'geneSetsDF' should not contain missing values.")
    }
    geneSetsDF$source <- as.character(geneSetsDF$source)
    numberOfSources <- dplyr::n_distinct(geneSetsDF$source)

    if(verbose) {
      message(
        paste0(
          numberOfSources,
          " different sources found."
        )
      )
    }

  }  else {
    geneSetsDF$source <- NA

    if(verbose) {
      message(
        paste0(
          "Assuming all gene-sets are from same source (multiple sources are supported)"
        )
      )
    }
  }

  ### Each gene set may appear once per source
  duplicatedKeys <- duplicated(geneSetsDF[, c('pathway', 'source')])
  if (any(duplicatedKeys)) {
    stop(
      "Each combination of 'pathway' and 'source' should appear only once in ",
      "'geneSetsDF', but these have several (different) rows: ",
      quoteAndCollapse(utils::head(unique(geneSetsDF$pathway[duplicatedKeys]), 5)),
      if (sum(duplicatedKeys) > 5) ", ..." else "", "."
    )
  }

  ### Gene sets to analyze
  if(!no_table){
    pathways <- unique(geneSetsDF$pathway)

    missingPathways <- setdiff(pathways, names(geneSetsList))
    if (length(missingPathways)) {
      stop(
        length(missingPathways), " pathway(s) of 'geneSetsDF' are not in ",
        "'geneSetsList', e.g. ",
        quoteAndCollapse(utils::head(missingPathways, 5)), "."
      )
    }
  } else{
    pathways <- names(geneSetsList)
  }

  ### Remove gene sets without genes
  emptySets <- pathways[lengths(geneSetsList[pathways]) == 0]
  if (length(emptySets)) {
    warning(
      length(emptySets), " gene set(s) contain no genes and were removed: ",
      quoteAndCollapse(utils::head(emptySets, 5)),
      if (length(emptySets) > 5) ", ..." else "", ".",
      call. = FALSE
    )
    pathways <- setdiff(pathways, emptySets)
    geneSetsDF <- geneSetsDF[geneSetsDF$pathway %in% pathways, , drop = FALSE]
  }

  if (!length(pathways)) {
    stop("No gene sets with genes are left to analyze.")
  }

  if(verbose) {
    message("Input is valid. Initializing gsList object...")
  }

  genesetsGenes <- geneSetsList[pathways]

  genesetsMatrix <-
    tibble::enframe(
      genesetsGenes,
      name = "geneset",
      value = "gene_ids"
    ) %>%
    tidyr::unnest(cols = c(gene_ids)) %>%
    dplyr::mutate(count = 1) %>%
    tidyr::pivot_wider(
      names_from = geneset,
      values_from = count,
      values_fill = list(count = 0)
    ) %>%
    tibble::column_to_rownames("gene_ids")

  res <- new("gsList", genesetsDataframe = genesetsMatrix)


  ### Drop the columns that would collide with those added downstream
  # Any other column is kept so it can be plotted via the 'plotFrom' argument
  # of e.g. 'plotEnrichments()' and 'plotHeatmap()'
  reservedColumns <- intersect(
    colnames(geneSetsDF),
    c('cluster', 'initial_cluster_label')
  )
  if( length(reservedColumns) ) {
    warning(paste0(
      "The following column(s) of 'geneSetsDF' are reserved for internal use ",
      'and were removed: ',
      paste0("'", reservedColumns, "'", collapse = ', '), '.'
    ))

    geneSetsDF <- geneSetsDF[
      , setdiff(colnames(geneSetsDF), reservedColumns), drop = FALSE
    ]
  }

  rownames(geneSetsDF) <- NULL
  res@enrichmentScore <- geneSetsDF

  numberOfSources <- dplyr::n_distinct(geneSetsDF$source)

  ### Record initialization metadata
  pkgVersion <- tryCatch(
    as.character(utils::packageVersion("geneSetSimplifyR")),
    error = function(e) NA_character_
  )
  res@metadata <- list(
    n_genes = nrow(res@genesetsDataframe),
    n_pathways = ncol(res@genesetsDataframe),
    n_sources = numberOfSources,
    package_version = pkgVersion,
    initialized_at = Sys.time()
  )

  if(verbose) {
    message(
      paste(
        "gsList object created with",
        dim(res@genesetsDataframe)[1],
        "genes and",
        dim(res@genesetsDataframe)[2],
        "gene sets",
        'from', numberOfSources, 'sources',
        sep = ' '
      )
    )
    message("Initialization of gsList object DONE")
  }

  return(res)
}
