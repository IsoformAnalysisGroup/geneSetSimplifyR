#' Example clustered and labeled gene sets
#'
#' A 'gsList' object made with `geneSetSimplifyR(exampleGeneSets,
#' exampleEnrichment, removeFirstWord = TRUE)`: 197 gene sets from 2 sources
#' ("Expression" and "Splicing") clustered at 10 resolutions. Use it to try the plotting and
#' table functions without running the workflow. The script that creates it is
#' in `inst/scripts/make-exampleGsList.R`.
#'
#' @format A 'gsList' object (see [gsList-class]).
#' @source Made from `exampleGeneSets` and `exampleEnrichment`.
"exampleGsList"

#' Example Enrichment Results
#'
#' fgsea results for the TCGA colon adenocarcinoma data analyzed in
#' Vitting-Seerup & Sandelin (2017), where gene-set enrichment was done both
#' for differentially expressed genes ("Expression") and differentially
#' spliced genes ("Splicing"). Only gene sets with an adjusted p-value below
#' 0.05 are included, a subset of 197 gene sets.
#'
#' @format A tibble with 206 rows (one per gene set and source) and 5 columns:
#' \describe{
#'   \item{source}{The analysis the gene set was found in: "Expression" or
#'     "Splicing".}
#'   \item{pathway}{The gene set name.}
#'   \item{enrichment_score}{Normalized enrichment score (NES) from fgsea.}
#'   \item{pval}{p-value from fgsea.}
#'   \item{padj}{Benjamini-Hochberg adjusted p-value from fgsea.}
#' }
#' @source Computed with `fgsea::fgseaMultilevel()` on MSigDB gene sets; see
#'   `inst/scripts/make-exampleEnrichment.R`. Vitting-Seerup K, Sandelin A.
#'   The Landscape of Isoform Switches in Human Cancers. Mol Cancer Res.
#'   2017;15(9):1206-1220.
"exampleEnrichment"

#' Example Gene Sets
#'
#' A named list of gene sets, where the names represent the gene set names and
#' each list element contains the corresponding genes for that gene set.
#'
#' @format A named list of gene sets, where the names represent the gene set
#'   names and each list element contains the corresponding genes for that
#'   gene set.
#' @source A subset of public data.
"exampleGeneSets"

#' Stop Words
#'
#' Words ignored when labeling clusters, e.g. "of", "and", "process". Bigrams
#' (two-word label candidates) containing one of these words are not used as
#' labels.
#'
#' @format A character vector.
#' @examples
#' head(stopWords)
#' # Add your own stop words
#' myStopWords <- c(stopWords, "regulation", "response")
"stopWords"
