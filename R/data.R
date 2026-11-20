#' Example clustered and labeled gene sets
#'
#' A 'gsList' object made with `geneSetSimplifyR(exampleGeneSets,
#' exampleEnrichment, removeFirstWord = TRUE)`: 196 gene sets from 2 sources
#' ("Expression" and "Splicing") clustered at 10 resolutions. Use it to try the
#' plotting and table functions without running the workflow. The script that
#' creates it is in `inst/scripts/make-exampleGsList.R`.
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
#' 0.05 are included, a subset of 196 gene sets.
#'
#' @format A tibble with 205 rows (one per gene set and source) and 5 columns:
#' \describe{
#'   \item{source}{The analysis the gene set was found in: "Expression" or
#'     "Splicing".}
#'   \item{pathway}{The gene set name.}
#'   \item{enrichment_score}{Normalized enrichment score (NES) from fgsea.}
#'   \item{pval}{p-value from fgsea.}
#'   \item{padj}{Benjamini-Hochberg adjusted p-value from fgsea.}
#' }
#' @source Computed with `fgsea::fgseaMultilevel()` and MSigDB v7.5.1 gene
#'   sets on the TCGA results of Vitting-Seerup K, Sandelin A. The Landscape
#'   of Isoform Switches in Human Cancers. Mol Cancer Res. 2017;15(9):1206-1220.
#'   \doi{10.1158/1541-7786.MCR-16-0459}. See
#'   `inst/scripts/make-exampleEnrichment.R`.
"exampleEnrichment"

#' Example Gene Sets
#'
#' The 196 gene sets that are significant in `exampleEnrichment`, as a named
#' list where the names are the gene-set names and each element holds the gene
#' symbols of that gene set.
#'
#' @format A named list of 196 character vectors of human gene symbols.
#' @source Human gene sets from the Molecular Signatures Database (MSigDB)
#'   v7.5.1 (R package msigdbr 7.5.1): the Hallmark collection, the C2 curated
#'   collection (chemical and genetic perturbations, Reactome and
#'   WikiPathways) and the C5 Gene Ontology collection. Each gene set is
#'   restricted to the genes tested in the enrichment analysis of
#'   `exampleEnrichment`. KEGG gene sets are left out because of their
#'   license. See `inst/scripts/make-exampleGeneSets.R`.
#'   <https://www.gsea-msigdb.org/gsea/msigdb>
"exampleGeneSets"

#' Stop Words
#'
#' Words ignored when labeling clusters, e.g. "of" and "and". Bigrams
#' (two-word label candidates) containing one of these words are not used as
#' labels.
#'
#' @format A character vector of 48 lower-case words: common English words
#'   (e.g. "of", "the"), gene-set collection prefixes (e.g. "gobp",
#'   "reactome"), the direction suffixes "up" and "dn" and the word "vs" from
#'   MSigDB gene-set names, author names that start MSigDB gene-set names, and
#'   other short fragments of gene-set names (e.g. "a5p4", "2").
#' @source Curated by the package authors.
#' @export
#' @examples
#' head(stopWords)
#' # Add your own stop words
#' myStopWords <- c(stopWords, "regulation", "response")
stopWords <- c(
    "gobp", "gomf", "hp", "gocc", "kegg", "gcnp", "reactome", "of", "to",
    "from", "by", "at", "in", "for", "with", "via", "into", "and", "or", "on",
    "up", "dn", "vs", "browne", "winter", "shepard", "baelde", "benporath",
    "blanco", "melo", "ginestier", "bertucci", "zhou", "mueller", "wp",
    "along", "a5p4", "amit", "d", "2", "turashvili", "rutella", "krige",
    "biocarta", "delys", "baris", "boylan", "the"
)
