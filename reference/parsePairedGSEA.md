# Parse pairedGSEA results into geneSetSimplifyR input

Converts `pairedGSEA::paired_ora()` output into the two inputs
[`geneSetSimplifyR()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/geneSetSimplifyR.md)
needs: `geneSetsList` and `geneSetsDF` (with `source`, `pathway`,
`enrichment_score`). Each of the three blocks (expression, splicing,
paired) is filtered by its own padj cutoff and added as a separate
`source`. Pathways passing multiple cutoffs appear once per source.

## Usage

``` r
parsePairedGSEA(
  oraResult,
  geneSets,
  cutoff_exp = 0.05,
  cutoff_splicing = 0.05,
  cutoff_paired = 0.05
)
```

## Arguments

- oraResult:

  Data frame from `pairedGSEA::paired_ora()`. Needs a `pathway` column
  plus `padj_<block>` and `enrichment_score_<block>` columns for each
  enabled block.

- geneSets:

  Named list of gene-set memberships (the same one passed to
  `pairedGSEA::paired_ora()`).

- cutoff_exp:

  padj cutoff for the expression block, or `NULL` to skip using
  expression. Default 0.05.

- cutoff_splicing:

  padj cutoff for the splicing block, or `NULL` to skip using splicing.
  Default 0.05.

- cutoff_paired:

  padj cutoff for the paired block, or `NULL` to skip it. Default 0.05.

## Value

A list with:

- `geneSetsList`:

  `geneSets` filtered to pathways passing at least one cutoff.

- `geneSetsDF`:

  Long-format data frame: `source`, `pathway`, `enrichment_score`,
  `padj` and `pval` (when `oraResult` has `pval_<block>` columns). One
  row per (pathway, source).

## Examples

``` r
# A small table in the format of pairedGSEA::paired_ora() output
pathways <- names(exampleGeneSets)[1:60]
oraResult <- data.frame(
  pathway = pathways,
  padj_expression = rep(c(0.01, 0.2), 30),
  enrichment_score_expression = seq(-2, 2, length.out = 60),
  padj_splicing = rep(c(0.2, 0.01, 0.01), 20),
  enrichment_score_splicing = seq(2, -2, length.out = 60),
  padj_paired = 1,
  enrichment_score_paired = 0
)
parsed <- parsePairedGSEA(
  oraResult       = oraResult,
  geneSets        = exampleGeneSets,
  cutoff_exp      = 0.05,
  cutoff_splicing = 0.05,
  cutoff_paired   = NULL
)
table(parsed$geneSetsDF$source)
#> 
#> expression   splicing 
#>         30         40 
gss <- geneSetSimplifyR(
  geneSetsList = parsed$geneSetsList,
  geneSetsDF   = parsed$geneSetsDF,
  resolution   = 1.1,
  verbose      = FALSE
)
```
