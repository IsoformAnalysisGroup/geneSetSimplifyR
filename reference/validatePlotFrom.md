# Validate the columns requested for plotting

Internal helper. Used by the plotting functions to check the `plotFrom`
argument before anything is summarized or drawn. A column is only usable
if it is present in the `enrichmentScore` slot, holds at least one
non-NA value, and is numeric.

## Usage

``` r
validatePlotFrom(geneSetsList, plotFrom)
```

## Arguments

- geneSetsList:

  A 'gsList' object.

- plotFrom:

  Name(s) of the column(s) in the `enrichmentScore` slot to plot, e.g.
  `'enrichment_score'`.

## Value

A single column name.

## Details

Only one column is plotted. Several can be supplied - the plot functions
declare their candidates as a default vector, the same way `summaryStat`
does - in which case the first usable one is returned. That way a
default of `c('nlog10pval', 'enrichment_score')` still works on a
`gsList` built without p-values, while a single explicitly named column
that cannot be plotted errors rather than silently showing something
else.
