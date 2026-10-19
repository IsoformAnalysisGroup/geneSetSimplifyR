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
#' required by some plotting functions (e.g. 'plotHeatmap()').
#' Optionally have a column named 'source' indicating the source of each gene set. This allows for
#' analysis of multi-conditional data throughout the entire workflow and is also supported by many
#' plotting functions. If a gene-set was found in multiple conditions there should be multiple rows
#' for the same 'pathway' indicating its respective source and corresponding enrichment_score.
#' See the vignette for examples.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @return 'gsList' object with a gene sets dataframe and enrichment scores for each gene set.
#' @export
#'
#' @examples
#' gsListObject <- initializeList(geneSetsList = exampleGeneSets, geneSetsDF = exampleEnrichment)
initializeList <- function(
    geneSetsList,
    geneSetsDF = NULL,
    verbose = TRUE
) {

  if(verbose) {
    message("Checking input...")
  }
  if (!is.list(geneSetsList)) {
    stop("'geneSetsList' should be a list")
  }

  if ( ! all(sapply(geneSetsList, is.character))) {
    stop("each element in the 'geneSetsList' should contain a character vector")
  }
  # Check input
  no_table <- FALSE
  if (is.null(geneSetsDF)) {
    geneSetsDF <- tibble::tibble(
      pathway = unique(names(geneSetsList)),
      enrichment_score = NA
    )
    no_table <- TRUE
  }

  # make unique
  geneSetsDF <- unique(geneSetsDF)

  # Check pathway
  if(! 'pathway' %in% colnames(geneSetsDF)) {
    stop("Column 'pathway' was not found in geneSetsDF.")
  }
  # Check enrichment_score
  if(! 'enrichment_score' %in% colnames(geneSetsDF)) {
    if( 'NES' %in% colnames(geneSetsDF) ) {
      message("enrichment_score was not found in the dataframe. Use NES instead.")
      geneSetsDF$enrichment_score <- geneSetsDF$NES
      geneSetsDF$NES <- NULL
    } else {
      geneSetsDF$enrichment_score <- NA
    }
  }

  # Check pvalues
  if(! 'pval' %in% colnames(geneSetsDF)) {
    if( 'pvalue' %in% colnames(geneSetsDF) ) {
      message("pval was not found in the dataframe. Use pvalue instead.")
      geneSetsDF$pval <- geneSetsDF$pvalue
      geneSetsDF$pvalue <- NULL
    } else {
      geneSetsDF$pval <- NA
    }
  }

  # Add the -log10 transformed p-values used by the plotting functions.
  # Done outside the block above so it also covers the (common) case where
  # 'pval' was supplied directly. p-values of 0 are floored at the smallest
  # representable double instead of becoming Inf; an absent p-value column
  # gives an all-NA column, which the plotting functions reject with an
  # informative error.
  if(! 'nlog10pval' %in% colnames(geneSetsDF)) {
    geneSetsDF$nlog10pval <- -log10(
      pmax(
        suppressWarnings(as.numeric(geneSetsDF$pval)),
        .Machine$double.xmin
      )
    )
  }

  if(!no_table){
    # Make pathway vector
    pathways <-
      geneSetsDF %>%
      dplyr::pull(pathway) %>%
      unique()

    if ( ! all(pathways %in% names(geneSetsList) ) ) {
      stop("not all genesets are in the provided database")
    }
  } else{
    pathways <- names(geneSetsList)
  }
  

  if(verbose) {
    message("Input is valid. Initializing gsList object...")
  }

  genesetsGenes <- geneSetsList[pathways]
  genesetsGenes <- lapply(genesetsGenes, unique)

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


  # Add enrichment score
  if( 'source' %in% colnames(geneSetsDF) ) {
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
