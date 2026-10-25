# Summarize enrichment per cluster and source

Internal helper. Summarizes the per-gene-set enrichment scores of
[`extractTable()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/extractTable.md)
into one value per cluster and source, using the statistic requested by
the plot function. Done here rather than reusing the `enrichment mean`
column of
[`extractTable()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/extractTable.md)
so the summary can follow `summaryStat`.

## Usage

``` r
summarizeEnrichment(
  geneSetsList,
  resolution = NULL,
  scoreColumn = "enrichment_score",
  summaryStat = "mean",
  removeClusterId = FALSE
)
```

## Arguments

- geneSetsList:

  A labeled 'gsList' object.

- resolution:

  `NULL` (the active resolution) or a clustered resolution.

- scoreColumn:

  Name of the column in the `enrichmentScore` slot to summarize, e.g.
  `'enrichment_score'` or `'nlog10pval'`. Should already have been
  checked with
  [`validatePlotFrom()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/validatePlotFrom.md)
  in the caller.

- summaryStat:

  One of 'mean', 'median' or 'max'. Must already be resolved to a single
  value (i.e. via [`match.arg()`](https://rdrr.io/r/base/match.arg.html)
  in the caller).

- removeClusterId:

  Logical. Drop the `"<id>: "` prefix from the labels.

## Value

A data.frame with columns `source`, `cluster`, `enrichmentSummary` and
`genesets` (the number of gene sets in the cluster and source, including
gene sets without a value), arranged by cluster.
