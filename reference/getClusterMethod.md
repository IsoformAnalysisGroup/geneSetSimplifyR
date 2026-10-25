# Resolve a clustering method specification

Turns the `algorithm` argument of
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
– a registered name, a `gsClusterMethod`, or a bare function – into a
`gsClusterMethod`.

## Usage

``` r
getClusterMethod(algorithm)
```

## Arguments

- algorithm:

  A registered method name, a `gsClusterMethod`, or a function following
  the cluster-method contract (see
  [`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md)).

## Value

A `gsClusterMethod`.

## Examples

``` r
getClusterMethod("louvain")
#> <gsClusterMethod> louvain
#>   Seurat::FindClusters, original Louvain (algorithm 1)
#>   requires: Seurat
```
