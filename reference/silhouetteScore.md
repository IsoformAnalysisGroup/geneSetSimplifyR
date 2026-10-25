# Mean silhouette width of a clustering

Internal helper for
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md).

## Usage

``` r
silhouetteScore(clusters, distMatrix)
```

## Arguments

- clusters:

  Integer cluster ids, one per gene set.

- distMatrix:

  Dense distance matrix between the gene sets.

## Value

The mean of the per-cluster average silhouette widths; 0 for a single
cluster and `NA` when the silhouette is not defined.
