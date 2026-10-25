# Apply a summary statistic while tolerating missing values

Internal helper. A column only has to be *partly* annotated to pass
[`validatePlotFrom()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/validatePlotFrom.md),
so NAs are dropped before summarizing rather than propagating to the
whole cluster. A cluster/source with no annotated gene set at all gives
`NA` instead of the `-Inf`/warning `max(numeric(0))` would give.

## Usage

``` r
summarizeScores(x, summaryFun)
```

## Arguments

- x:

  Numeric vector.

- summaryFun:

  Function to apply, e.g. `mean`.

## Value

A single number, or `NA_real_` if `x` holds no non-NA value.
