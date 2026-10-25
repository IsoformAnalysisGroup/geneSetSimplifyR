# Register a clustering method for use by name

Adds a clustering backend to the package registry so it can be selected
with `clusterGeneSets(algorithm = "<name>")`. Registering is optional –
a `gsClusterMethod` (or a bare function) can also be passed straight to
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md).

## Usage

``` r
registerClusterMethod(
  name,
  fn,
  description = NA_character_,
  packages = character(0),
  usesResolution = TRUE,
  overwrite = FALSE
)
```

## Arguments

- name:

  Name to register the method under.

- fn:

  Either a `gsClusterMethod` (from
  [`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md))
  or a function following the same contract; a bare function is wrapped
  automatically.

- description:

  One-line description; ignored when `fn` is already a
  `gsClusterMethod`.

- packages:

  Packages the method requires; ignored when `fn` is already a
  `gsClusterMethod`.

- usesResolution:

  Whether the method's output depends on `resolution`; ignored when `fn`
  is already a `gsClusterMethod`. See
  [`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md).

- overwrite:

  Whether to replace an already-registered method of the same name.
  Defaults to `FALSE`.

## Value

The registered `gsClusterMethod`, invisibly.

## See also

[`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md),
[`listClusterMethods()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/listClusterMethods.md)

## Examples

``` r
registerClusterMethod(
  "leiden_cpm",
  function(snn, resolution, verbose = FALSE, ...) {
    g <- snnToIgraph(snn)
    igraph::membership(igraph::cluster_leiden(
      g,
      objective_function = "CPM",
      resolution = resolution,
      weights = igraph::E(g)$weight
    ))
  },
  description = "igraph Leiden with the CPM objective",
  packages = "igraph",
  overwrite = TRUE
)
listClusterMethods()
#>                  name
#> 1           apcluster
#> 2          binary_cut
#> 3       igraph_leiden
#> 4      igraph_louvain
#> 5              kmeans
#> 6              leiden
#> 7          leiden_cpm
#> 8             louvain
#> 9  louvain.multilevel
#> 10                slm
#>                                                                                               description
#> 1                                Affinity propagation (apcluster::apcluster) on the SNN similarity matrix
#> 2                                             simplifyEnrichment::binary_cut on the SNN similarity matrix
#> 3                                          igraph::cluster_leiden on the SNN graph (modularity objective)
#> 4                                                                igraph::cluster_louvain on the SNN graph
#> 5  k-means on SNN similarity profiles; k chosen by elbow/knee on the WSS curve (as in simplifyEnrichment)
#> 6                                                Seurat::FindClusters, Leiden via leidenalg (algorithm 4)
#> 7                                                                    igraph Leiden with the CPM objective
#> 8                                                    Seurat::FindClusters, original Louvain (algorithm 1)
#> 9                                  Seurat::FindClusters, Louvain with multilevel refinement (algorithm 2)
#> 10                                                 Seurat::FindClusters, smart local moving (algorithm 3)
#>              packages uses_resolution
#> 1           apcluster           FALSE
#> 2  simplifyEnrichment           FALSE
#> 3              igraph            TRUE
#> 4              igraph            TRUE
#> 5                               FALSE
#> 6      Seurat, leiden            TRUE
#> 7              igraph            TRUE
#> 8              Seurat            TRUE
#> 9              Seurat            TRUE
#> 10             Seurat            TRUE

geneSetsList <- initializeList(exampleGeneSets, exampleEnrichment,
                               verbose = FALSE)
gsListObject <- clusterGeneSets(geneSetsList, algorithm = "leiden_cpm",
                                resolution = 0.05, verbose = FALSE)
```
