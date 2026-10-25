# Normalize cluster factor level order in a clustered gsList

Rebuilds every `GSEA_snn_res.<res>` metadata column so its factor levels
are in numeric order (0, 1, 2, ...), starting at 0. Newer objects get
this automatically inside
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md);
use this to migrate objects clustered before that fix without
re-clustering. Cluster assignments are unchanged – only the level order
is corrected.

## Usage

``` r
normalizeClusterLevels(geneSetsList)
```

## Arguments

- geneSetsList:

  A clustered 'gsList' object.

## Value

The 'gsList' with normalized cluster factor level order.

## Examples

``` r
gsListObject <- normalizeClusterLevels(exampleGsList)
levels(gsListObject@cluster$GSEA_snn_res.1.1)
#> [1] "0" "1" "2" "3" "4" "5"
```
