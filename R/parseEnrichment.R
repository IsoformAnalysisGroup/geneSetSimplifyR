#' Parse pairedGSEA results into geneSetSimplifyR input
#'
#' @description Converts `pairedGSEA::paired_ora()` output into the two inputs
#'   `geneSetSimplifyR()` needs: `geneSetsList` and `geneSetsDF` (with `source`,
#'   `pathway`, `enrichment_score`). Each of the three blocks (expression,
#'   splicing, paired) is filtered by its own padj cutoff and added as a
#'   separate `source`. Pathways passing multiple cutoffs appear once per source.
#'
#' @param oraResult Data frame from `pairedGSEA::paired_ora()`. Needs a
#'   `pathway` column plus `padj_<block>` and `enrichment_score_<block>`
#'   columns for each enabled block.
#' @param geneSets Named list of gene-set memberships (the same one passed to
#'   `pairedGSEA::paired_ora()`).
#' @param cutoff_exp padj cutoff for the expression block, or `NULL` to skip
#'   using expression. Default 0.05.
#' @param cutoff_splicing padj cutoff for the splicing block, or `NULL` to
#'   skip using splicing. Default 0.05.
#' @param cutoff_paired padj cutoff for the paired block, or `NULL` to skip
#'   it. Default 0.05.
#'
#' @return A list with:
#'   \describe{
#'     \item{`geneSetsList`}{`geneSets` filtered to pathways passing at least
#'       one cutoff.}
#'     \item{`geneSetsDF`}{Long-format data frame: `source`, `pathway`,
#'       `enrichment_score`. One row per (pathway, source).}
#'   }
#'
#' @examples
#' \dontrun{
#' parsed <- parsePairedGSEA(
#'   oraResult       = oraResult,
#'   geneSets        = mSigDB,
#'   cutoff_exp      = 0.05,
#'   cutoff_splicing = 0.05,
#'   cutoff_paired   = 0.05
#' )
#' gss <- geneSetSimplifyR(
#'   geneSetsList = parsed$geneSetsList,
#'   geneSetsDF   = parsed$geneSetsDF
#' )
#' }
#'
#' @export
parsePairedGSEA <- function(
    oraResult,
    geneSets,
    cutoff_exp      = 0.05,
    cutoff_splicing = 0.05,
    cutoff_paired   = 0.05
) {
  if (missing(oraResult) || is.null(oraResult)) {
    stop("'oraResult' is required.")
  }
  if (missing(geneSets) || is.null(geneSets)) {
    stop("'geneSets' is required -- pass the same gene_sets list you used when running pairedGSEA::paired_ora().")
  }
  if (!is.list(geneSets) || is.null(names(geneSets))) {
    stop("'geneSets' should be a named list (pathway name -> character vector of gene IDs).")
  }

  oraDF <- as.data.frame(oraResult, stringsAsFactors = FALSE)

  if (!"pathway" %in% colnames(oraDF)) {
    stop("'oraResult' must contain a 'pathway' column.")
  }

  blocks <- list(
    list(cutoff = cutoff_exp,      tag = "expression", source = "expression"),
    list(cutoff = cutoff_splicing, tag = "splicing",   source = "splicing"),
    list(cutoff = cutoff_paired,   tag = "paired",     source = "paired")
  )

  long_rows <- list()
  for (b in blocks) {
    if (is.null(b$cutoff)) next
    if (!is.numeric(b$cutoff) || length(b$cutoff) != 1 || b$cutoff <= 0 || b$cutoff > 1) {
      stop(sprintf("Cutoff for '%s' must be a single number in (0, 1] (got: %s).",
                   b$source, paste(b$cutoff, collapse = ", ")))
    }

    padj_col <- paste0("padj_", b$tag)
    es_col   <- paste0("enrichment_score_", b$tag)

    if (!padj_col %in% colnames(oraDF)) {
      stop(sprintf(
        "Column '%s' not found in oraResult -- cannot apply cutoff for '%s'. Set cutoff_%s = NULL to skip this block.",
        padj_col, b$source, b$tag
      ))
    }
    if (!es_col %in% colnames(oraDF)) {
      stop(sprintf("Column '%s' not found in oraResult.", es_col))
    }

    keep <- !is.na(oraDF[[padj_col]]) & oraDF[[padj_col]] < b$cutoff
    if (any(keep)) {
      long_rows[[b$source]] <- data.frame(
        source           = b$source,
        pathway          = oraDF$pathway[keep],
        enrichment_score = oraDF[[es_col]][keep],
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(long_rows)) {
    stop("No pathways passed any of the enabled cutoffs. Try relaxing the cutoffs or check the input.")
  }

  geneSetsDF <- do.call(rbind, long_rows)
  rownames(geneSetsDF) <- NULL

  # Drop pathways not in geneSets (warn the user once with a sample)
  uniquePaths   <- unique(geneSetsDF$pathway)
  missingPaths  <- setdiff(uniquePaths, names(geneSets))
  if (length(missingPaths)) {
    warning(sprintf(
      "%d pathway(s) in oraResult are not in 'geneSets' and will be dropped (first 5: %s).",
      length(missingPaths),
      paste(utils::head(missingPaths, 5), collapse = ", ")
    ), call. = FALSE)
    geneSetsDF  <- geneSetsDF[geneSetsDF$pathway %in% names(geneSets), , drop = FALSE]
    uniquePaths <- unique(geneSetsDF$pathway)
  }

  if (!length(uniquePaths)) {
    stop("After dropping pathways not in 'geneSets', no pathways remain.")
  }

  list(
    geneSetsList = geneSets[uniquePaths],
    geneSetsDF   = geneSetsDF
  )
}


#' Parse fgsea results into geneSetSimplifyR input
#'
#' @description Converts [fgsea::fgsea()] output into the two inputs
#'   `geneSetSimplifyR()` needs: `geneSetsList` and `geneSetsDF` (with
#'   `source`, `pathway`, `enrichment_score`). To combine several fgsea runs
#'   into one multi-source analysis, call this function once per run with a
#'   distinct `source` label and `rbind()` the resulting `geneSetsDF`s.
#'
#' @param fgseaResult Output of `fgsea::fgsea()`. Must contain `pathway`,
#'   `padj`, and `NES` columns.
#' @param geneSets Named list of gene-set memberships (the same one passed to
#'   `fgsea::fgsea()`).
#' @param cutoff padj threshold for filtering pathways. Default 0.05.
#' @param source Label written to the `source` column of `geneSetsDF`. Use a
#'   distinct value per run when combining multiple analyses. Default
#'   `"fgsea"`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`geneSetsList`}{`geneSets` filtered to pathways passing the
#'       cutoff.}
#'     \item{`geneSetsDF`}{Long-format data frame: `source`, `pathway`,
#'       `enrichment_score` (fgsea's NES). One row per pathway.}
#'   }
#'
#' @examples
#' \dontrun{
#' parsed <- parseFgsea(
#'   fgseaResult = fgseaRes,
#'   geneSets    = mSigDB,
#'   cutoff      = 0.05,
#'   source      = "treatmentA"
#' )
#' gss <- geneSetSimplifyR(
#'   geneSetsList = parsed$geneSetsList,
#'   geneSetsDF   = parsed$geneSetsDF
#' )
#' }
#'
#' @export
parseFgsea <- function(
    fgseaResult,
    geneSets,
    cutoff = 0.05,
    source = "fgsea"
) {
  if (missing(fgseaResult) || is.null(fgseaResult)) {
    stop("'fgseaResult' is required.")
  }
  if (missing(geneSets) || is.null(geneSets)) {
    stop("'geneSets' is required -- pass the same pathways list you used when running fgsea::fgsea().")
  }
  if (!is.list(geneSets) || is.null(names(geneSets))) {
    stop("'geneSets' should be a named list (pathway name -> character vector of gene IDs).")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1 || cutoff <= 0 || cutoff > 1) {
    stop(sprintf("'cutoff' must be a single number in (0, 1] (got: %s).",
                 paste(cutoff, collapse = ", ")))
  }
  if (!is.character(source) || length(source) != 1 || !nzchar(source)) {
    stop("'source' must be a single non-empty character string.")
  }

  df <- as.data.frame(fgseaResult, stringsAsFactors = FALSE)

  requiredCols <- c("pathway", "padj", "NES")
  missingCols  <- setdiff(requiredCols, colnames(df))
  if (length(missingCols)) {
    if (length(missingCols) == 1 & 
          missingCols == "NES" & 
          "foldEnrichment" %in% colnames(df)){
      df$NES <- df$foldEnrichment
    } else {
      stop(sprintf("'fgseaResult' is missing required column(s): %s.",
           paste(missingCols, collapse = ", ")))
    }
  }

  keep <- !is.na(df$padj) & df$padj < cutoff
  if (!any(keep)) {
    stop(sprintf("No pathways passed the cutoff (padj < %s).", cutoff))
  }

  geneSetsDF <- data.frame(
    source           = source,
    pathway          = df$pathway[keep],
    enrichment_score = df$NES[keep],
    stringsAsFactors = FALSE
  )

  # Drop pathways not in geneSets (warn once with a sample)
  uniquePaths  <- unique(geneSetsDF$pathway)
  missingPaths <- setdiff(uniquePaths, names(geneSets))
  if (length(missingPaths)) {
    warning(sprintf(
      "%d pathway(s) in fgseaResult are not in 'geneSets' and will be dropped (first 5: %s).",
      length(missingPaths),
      paste(utils::head(missingPaths, 5), collapse = ", ")
    ), call. = FALSE)
    geneSetsDF  <- geneSetsDF[geneSetsDF$pathway %in% names(geneSets), , drop = FALSE]
    uniquePaths <- unique(geneSetsDF$pathway)
  }

  if (!length(uniquePaths)) {
    stop("After dropping pathways not in 'geneSets', no pathways remain.")
  }

  list(
    geneSetsList = geneSets[uniquePaths],
    geneSetsDF   = geneSetsDF
  )
}
