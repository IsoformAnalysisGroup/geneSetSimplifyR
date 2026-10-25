# Plot enrichment as heatmap

Creates a heatmap summarizing the enrichment of each gene-set cluster in
each source (useful for datasets with many data sources).

## Usage

``` r
plotHeatmap(
  geneSetsList,
  resolution = NULL,
  plotFrom = c("nlog10pval", "enrichment_score"),
  summaryStat = c("median", "max", "mean"),
  removeClusterId = FALSE,
  highlightClusterNo = NULL,
  transpose = FALSE,
  cutree_rows = NA,
  cutree_cols = NA,
  drawPlot = TRUE,
  scaleCols = "auto"
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

- plotFrom:

  Name of the column of the 'enrichmentScore' slot to visualize, e.g.
  'nlog10pval' or 'enrichment_score'. Only numeric columns that are
  annotated (i.e. not all NA) can be plotted - which columns are
  available depends on what was supplied to 'geneSetsDF' in
  'initializeList()'. Only one column is plotted; if several are
  supplied the first plottable one is used (the default therefore falls
  back to 'enrichment_score' when no p-values were supplied).

- summaryStat:

  How the values of the gene sets within a cluster are summarized into
  the value shown in each heatmap cell. One of 'median' (default), 'max'
  or 'mean'. Note that 'max' is the largest signed value, so for signed
  scores (e.g. NES) a strongly negative cluster will be summarized by
  its least negative gene set.

- removeClusterId:

  Logical, if TRUE the cluster number will be removed from each cluster
  label. Default: FALSE.

- highlightClusterNo:

  A vector with the cluster numbers of the clusters (at the resolution)
  to highlight (in bold) in the plot.

- transpose:

  A logical indicating whether to transpose the heatmap (switch rows and
  cols)

- cutree_rows:

  pheatmap control parameter. Number of clusters the rows are divided
  into

- cutree_cols:

  pheatmap control parameter. Number of clusters the columns are divided
  into

- drawPlot:

  Whether to draw the plot. Default to TRUE.

- scaleCols:

  Whether to z-score each cluster across the sources before plotting.
  One of `"auto"` (default), `TRUE` or `FALSE`. `"auto"` scales only
  when there are at least three sources: with one source there is
  nothing to compare, and with two sources the z-scores can only be
  -0.71, 0 or 0.71. Scaling is always across sources, also when
  `transpose = TRUE`.

## Value

The pheatmap object, invisibly. Use `drawPlot = FALSE` to get it without
drawing, e.g. to arrange it with other plots.

## Details

Cells of clusters that have no gene set from a source are shown in gray.
Rows and columns are ordered by hierarchical clustering when there are
at least two sources (absent cells count as 0 for the clustering).

## Examples

``` r
plotHeatmap(exampleGsList)

plotHeatmap(exampleGsList, plotFrom = 'enrichment_score', summaryStat = 'mean')
```
