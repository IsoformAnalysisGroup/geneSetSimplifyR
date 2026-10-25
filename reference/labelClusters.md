# Label clusters in clustered gene sets list object

Label clusters in clustered gene sets list object

## Usage

``` r
labelClusters(
  geneSetsList,
  stopWords = geneSetSimplifyR::stopWords,
  nLabels = 10,
  removeFirstWord = FALSE,
  verbose = TRUE
)
```

## Arguments

- geneSetsList:

  'gsList' object that has been clustered with the 'clusterGeneSets()'
  function.

- stopWords:

  Vector with stop words. The default vector can be viewed with
  data(stopWords). To add new stop words 'stopWords \<-
  c(geneSetSimplifyR::stopWords, 'new word 1', 'new word 2', ...)

- nLabels:

  Number of top labels to keep. Default to 10

- removeFirstWord:

  Logical. If TRUE, strip the first word (everything before the first
  underscore) from each gene-set name before generating bigram tokens.
  Useful for MSigDB-style names where the first token is an author or
  database prefix (e.g. "GOBP", "GOMF", "HALLMARK"). Default: FALSE.

- verbose:

  a logic indicating whether to print progress statements. Defaults to
  TRUE.

## Value

'gsList' object with labeled clusters for all clustered resolutions.

## Examples

``` r
# Relabel the shipped example, e.g. with extra stop words
gsListObject <- labelClusters(
  geneSetsList = exampleGsList,
  stopWords = c(stopWords, "regulation"),
  verbose = FALSE
)
```
