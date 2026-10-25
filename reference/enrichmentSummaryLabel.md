# Axis/legend text describing an enrichment summary

Internal helper.

## Usage

``` r
enrichmentSummaryLabel(summaryStat = "mean", scoreColumn = "enrichment_score")
```

## Arguments

- summaryStat:

  One of 'mean', 'median' or 'max'.

- scoreColumn:

  Name of the summarized column. `'enrichment_score'` is spelled out as
  "Enrichment"; any other column is named as-is so the axis says what is
  actually plotted.

## Value

A single string.
