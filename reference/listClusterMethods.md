# List the available clustering methods

Returns the clustering backends
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
can be pointed at by name, including any added with
[`registerClusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/registerClusterMethod.md).

## Usage

``` r
listClusterMethods()
```

## Value

A `data.frame` with columns `name`, `description`, `packages` and
`uses_resolution` (`FALSE` for methods that choose the number of
clusters themselves, which
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
therefore runs only once).

## See also

[`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md),
[`registerClusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/registerClusterMethod.md)

## Examples

``` r
listClusterMethods()
#>                 name
#> 1          apcluster
#> 2         binary_cut
#> 3      igraph_leiden
#> 4     igraph_louvain
#> 5             kmeans
#> 6             leiden
#> 7            louvain
#> 8 louvain.multilevel
#> 9                slm
#>                                                                                              description
#> 1                               Affinity propagation (apcluster::apcluster) on the SNN similarity matrix
#> 2                                            simplifyEnrichment::binary_cut on the SNN similarity matrix
#> 3                                         igraph::cluster_leiden on the SNN graph (modularity objective)
#> 4                                                               igraph::cluster_louvain on the SNN graph
#> 5 k-means on SNN similarity profiles; k chosen by elbow/knee on the WSS curve (as in simplifyEnrichment)
#> 6                                               Seurat::FindClusters, Leiden via leidenalg (algorithm 4)
#> 7                                                   Seurat::FindClusters, original Louvain (algorithm 1)
#> 8                                 Seurat::FindClusters, Louvain with multilevel refinement (algorithm 2)
#> 9                                                 Seurat::FindClusters, smart local moving (algorithm 3)
#>             packages uses_resolution
#> 1          apcluster           FALSE
#> 2 simplifyEnrichment           FALSE
#> 3             igraph            TRUE
#> 4             igraph            TRUE
#> 5                              FALSE
#> 6     Seurat, leiden            TRUE
#> 7             Seurat            TRUE
#> 8             Seurat            TRUE
#> 9             Seurat            TRUE
```
