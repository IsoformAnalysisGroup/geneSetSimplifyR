# Cluster names to highlight in a plot

Internal helper. Translates cluster numbers into the cluster names shown
in the plot of the given resolution.

## Usage

``` r
highlightedClusterNames(
  geneSetsList,
  resolution,
  highlightClusterNo,
  removeClusterId
)
```

## Arguments

- geneSetsList:

  A labeled 'gsList' object.

- resolution:

  `NULL` (the active resolution) or a clustered resolution.

- highlightClusterNo:

  Numeric vector of cluster numbers, or `NULL`.

- removeClusterId:

  Logical. Whether the plot drops the `"<id>: "` prefix.

## Value

Character vector of cluster names (empty when nothing to highlight).
