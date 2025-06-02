#' Initialize gene sets list
#'
#' @description Initialize a 'gsList' object from a dataframe with gene sets and enrichment scores and a gene sets database.
#'
#' @param geneSetsList A list containing the gene sets to be analyzed.
#' The name of each element in the list should be a gene set and each element should contain
#' the genes belonging to that gene set.
#' @param geneSetsDF Optional data.frame with gene sets and enrichment scores that can be usefull
#' for downstream visualization (both for effect size and direction and if multiple conditions are analyzed).
#' Should have a column named 'pathway' containing the gene sets and a column
#' named 'enrichment_score' or 'NES' with the enrichment of each gene set.
#' Optinally have a column named 'source' indicating the source of each gene set. This allows for
#' analysis of multi-condtional data throughput the entire workflow and is also supported by many
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
  if (is.null(geneSetsDF)) {
    geneSetsDF <- tibble::tibble(
      pathway = unique(names(geneSetsList)),
      enrichment_score = NA
    )
  }

  # make unique
  geneSetsDF <- unique(geneSetsDF)
  if(! 'enrichment_score' %in% colnames(geneSetsDF)) {
    if( 'NES' %in% colnames(geneSetsDF) ) {
      geneSetsDF$enrichment_score <- geneSetsDF$NES
      geneSetsDF$NES <- NULL
    } else {
      geneSetsDF$enrichment_score <- NA
    }
  }

  # Make pathway vector
  pathways <-
    geneSetsDF %>%
    dplyr::pull(pathway) %>%
    unique()

  if ( ! all(pathways %in% names(geneSetsList) ) ) {
    stop("not all genesets are in the provided database")
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

  ### Subset (to ensure no confounding columns downstream)
  geneSetsDF <- geneSetsDF[,c('source','pathway','enrichment_score')]

  res@enrichmentScore <- geneSetsDF

  numberOfSources <- dplyr::n_distinct(geneSetsDF$source)

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
