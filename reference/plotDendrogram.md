# Plot dendrogram of cluster averages

A plot showing the hierarchical clustering of the cluster averages
(dendrogram). A tool made for showing how clusters (not individual
gene-sets) relate to each other. The distance between two clusters is
the average SNN distance (1 - SNN similarity) between their gene sets,
and clusters are joined by complete linkage.

## Usage

``` r
plotDendrogram(
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

A ggplot object. The underlying `hclust` object is stored in its
`hclust` attribute.

## Examples

``` r
plotDendrogram(exampleGsList)

hc <- attr(plotDendrogram(exampleGsList), "hclust")
```
