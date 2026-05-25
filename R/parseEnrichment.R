#' Parse pairedGSEA over-representation results into geneSetSimplifyR input
#'
#' @description Converts the output of `pairedGSEA::paired_ora()` into the two
#'   inputs `geneSetSimplifyR()` expects: a `geneSetsList` (named list of gene
#'   memberships) and a long-format `geneSetsDF` with `source`, `pathway`,
#'   `enrichment_score`. Each of the three pairedGSEA result blocks
#'   (expression, splicing, paired) is filtered independently by its own padj
#'   cutoff and contributes its surviving pathways to `geneSetsDF` as a
#'   separate `source`. A pathway that passes multiple cutoffs appears in
#'   `geneSetsDF` multiple times (once per source) -- the standard long
#'   format used elsewhere in the package.
#'
#' @param oraResult The data frame (or `DFrame`) returned by
#'   `pairedGSEA::paired_ora()`. Must contain a `pathway` column and, for each
#'   enabled block, the `padj_<block>` and `enrichment_score_<block>` columns.
#' @param geneSets Named list of gene-set memberships -- the same `gene_sets`
#'   used when running `pairedGSEA::paired_ora()` (e.g. from
#'   `pairedGSEA::prepare_msigdb()` or [getGeneSets()]).
#' @param cutoff_exp Numeric in (0, 1] or `NULL`. padj cutoff for the
#'   `_expression` block. `NULL` disables this block. Default 0.05.
#' @param cutoff_splicing Numeric in (0, 1] or `NULL`. padj cutoff for the
#'   `_splicing` block. `NULL` disables this block. Default 0.05.
#' @param cutoff_paired Numeric in (0, 1] or `NULL`. padj cutoff for the
#'   `_paired` block. `NULL` disables this block. Default 0.05.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`geneSetsList`}{Subset of `geneSets` restricted to pathways that
#'       survived at least one of the enabled cutoffs.}
#'     \item{`geneSetsDF`}{Long-format data frame with columns `source`
#'       (\code{"expression"} / \code{"splicing"} / \code{"paired"}),
#'       `pathway`, `enrichment_score`. One row per (pathway, source) pair.}
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
