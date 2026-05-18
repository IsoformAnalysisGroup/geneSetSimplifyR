#' Parse fgsea enrichment results (work-in-progress)
#'
#' Reads an fgsea result table and returns a data frame with `pathway`
#' and `NES` columns for gene sets passing the FDR cutoff.
#'
#' @param input data.frame of fgsea output (must contain `padj`, `pathway`, `NES`).
#' @param geneList Optional gene list (currently unused; reserved for future use).
#' @param FDRCutoff Numeric. FDR threshold for filtering pathways.
#'
#' @return A data frame with columns `pathway` and `NES`.
#'
#' @keywords internal
#' @noRd
read_fgsea <- function(input, geneList, FDRCutoff = 0.05) {
  input %>%
    dplyr::filter(padj <= FDRCutoff) %>%
    dplyr::select(pathway, NES)
}
