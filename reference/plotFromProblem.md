# Why a column cannot be plotted

Internal helper for
[`validatePlotFrom()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/validatePlotFrom.md).

## Usage

``` r
plotFromProblem(enrichmentDF, column)
```

## Arguments

- enrichmentDF:

  The `enrichmentScore` slot of a 'gsList' object.

- column:

  Name of a single column.

## Value

`NA_character_` if the column can be plotted, otherwise a one-line
description of the problem.

## Details

The three checks are done in this order on purpose: a score that was
never supplied to
[`initializeList()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/initializeList.md)
is stored as a logical `NA` column, so checking for "all NA" before "is
numeric" gives the informative reason rather than a confusing complaint
about the column type.
