#' Example clustered and labeled gene sets
#'
#' A 'gsList' object made with `geneSetSimplifyR(exampleGeneSets,
#' exampleEnrichment)`: 197 gene sets from 2 sources ("Expression" and
#' "Splicing") clustered at 10 resolutions. Use it to try the plotting and
#' table functions without running the workflow. The script that creates it is
#' in `inst/scripts/make-exampleGsList.R`.
#'
#' @format A 'gsList' object (see [gsList-class]).
#' @source Made from `exampleGeneSets` and `exampleEnrichment`.
"exampleGsList"

#' Example Enrichment Results
#'
#' A subset of the colorectal data where gene-set enrichment was done both for differentially expressed genes ("Expression") and differentially spliced genes ("Splicing").
#'
#' @format A data frame with 3 columns, source, pathway, and enrichment_score.
#' \describe{
#'   \item{source}{The source/group used when applying GSEA.}
#'   \item{pathway}{The pathway name when applying GSEA.}
#'   \item{enrichment_score}{Enrichment scores from GSEA results.}
#' }
#' @source A subset of public data.
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
