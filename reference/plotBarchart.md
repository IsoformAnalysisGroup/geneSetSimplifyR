# Plot clusters as a bar chart

Creates a bar chart showing the number of gene sets in each cluster, per
source when there are several.

## Usage

``` r
plotBarchart(
  geneSetsList,
  resolution = NULL,
  removeClusterId = FALSE,
  highlightClusterNo = NULL
)
```

## Arguments

- geneSetsList:

  'gsList' object with clusters and labels.

- resolution:

  The cluster resolution to plot. Default `NULL` uses the active
  resolution: the optimal resolution found by
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md),
  or the one last set with
  [`updateResolution()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateResolution.md)
  or
  [`updateClusterLabels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateClusterLabels.md).

- removeClusterId:

  Logical, if TRUE the cluster number will be removed from each cluster
  label. Default: FALSE.

- highlightClusterNo:

  A vector with the cluster numbers of the clusters (at the resolution)
  to highlight (in bold) in the plot.

## Value

A ggplot object.

## Examples

``` r
plotBarchart(exampleGsList)

plotBarchart(exampleGsList, resolution = 0.6, highlightClusterNo = 1)
```
