## Creates data/exampleGsList.rda, the clustered and labeled example object.
##
## Run from the root of the package source with the package installed (or
## loaded with pkgload::load_all()):
##   Rscript inst/scripts/make-exampleGsList.R
##
## Input: the shipped example data 'exampleGeneSets' (gene sets) and
## 'exampleEnrichment' (enrichment scores of the gene sets in two sources,
## "Expression" and "Splicing", a subset of a colorectal cancer data set).

if (!"geneSetSimplifyR" %in% loadedNamespaces()) {
  library(geneSetSimplifyR)
}

set.seed(42)
exampleGsList <- geneSetSimplifyR(
  geneSetsList = exampleGeneSets,
  geneSetsDF = exampleEnrichment,
  verbose = FALSE,
  # The gene-set names start with the name of their author, e.g.
  # SCHLOSSER_SERUM_RESPONSE_DN; leave it out of the labels
  removeFirstWord = TRUE
)

## Keep the shipped object small: the LSI feature loadings and the raw SVD
## output are not used by any function of the package.
exampleGsList@cluster@reductions$lsi@feature.loadings <- matrix(numeric(0), 0, 0)
exampleGsList@cluster@reductions$lsi@misc <- list()

save(
  exampleGsList,
  file = file.path("data", "exampleGsList.rda"),
  compress = "xz"
)
