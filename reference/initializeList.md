# Initialize gene sets list

Initialize a 'gsList' object from a dataframe with gene sets and
enrichment scores and a gene sets database.

## Usage

``` r
initializeList(geneSetsList, geneSetsDF = NULL, verbose = TRUE)
```

## Arguments

- geneSetsList:

  A list containing the gene sets to be analyzed. The name of each
  element in the list should be a gene set and each element should
  contain the genes belonging to that gene set. Alternatively, the path
  to a GMT file with the gene sets (read with
  [`fgsea::gmtPathways()`](https://rdrr.io/pkg/fgsea/man/gmtPathways.html)).

- geneSetsDF:

  Optional data.frame with gene sets and enrichment scores that can be
  useful for downstream visualization (both for effect size and
  direction and if multiple conditions are analyzed). Should have a
  column named 'pathway' containing the gene sets and a column named
  'enrichment_score' or 'NES' with the enrichment of each gene set.
  Optionally have a column named 'pval' or 'pvalue' with the p-value of
  each gene set; a 'nlog10pval' column with the -log10 transformed
  p-values is added automatically and is used by default by the plotting
  functions that summarize scores (e.g. 'plotHeatmap()'). p-values must
  lie between 0 and 1. A p-value of 0 is transformed as the smallest
  positive double (giving about 307.7) rather than infinity. Optionally
  have a column named 'source' indicating the source of each gene set.
  This allows for analysis of multi-conditional data throughout the
  entire workflow and is also supported by many plotting functions. If a
  gene-set was found in multiple conditions there should be multiple
  rows for the same 'pathway' indicating its respective source and
  corresponding enrichment_score. Each combination of 'pathway' and
  'source' may appear only once. Any other numeric column is kept and
  can be plotted via the 'plotFrom' argument of the plotting functions.
  List columns (e.g. 'leadingEdge' from fgsea) are dropped. See the
  vignette for examples.

- verbose:

  a logic indicating whether to print progress statements. Defaults to
  TRUE.

## Value

'gsList' object with a gene sets dataframe and enrichment scores for
each gene set.

## Details

Gene sets without any gene are removed with a warning. Missing gene ids
(`NA` or empty strings) are ignored.

## Examples

``` r
gsListObject <- initializeList(geneSetsList = exampleGeneSets, geneSetsDF = exampleEnrichment)
#> Checking input...
#> 2 different sources found.
#> Input is valid. Initializing gsList object...
#> gsList object created with 11210 genes and 197 gene sets from 2 sources
#> Initialization of gsList object DONE
gsListObject
#> A gsList object with 197 gene sets from 2 source(s).
#> Metadata:
#>   n pathways: 197

# Gene sets can also be read from a GMT file
gmtFile <- tempfile(fileext = ".gmt")
writeLines(
  vapply(names(exampleGeneSets), function(name) {
    paste(c(name, "", exampleGeneSets[[name]]), collapse = "\t")
  }, character(1)),
  gmtFile
)
gsListObject <- initializeList(gmtFile, exampleEnrichment, verbose = FALSE)
```
