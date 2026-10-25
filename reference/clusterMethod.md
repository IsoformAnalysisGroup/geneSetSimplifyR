# Define a clustering method for `clusterGeneSets()`

Wraps a community-detection function into the object
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
expects. The function receives the shared-nearest- neighbour (SNN) graph
built by
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
and a single resolution, and returns one cluster id per gene set.

## Usage

``` r
clusterMethod(
  fn,
  name,
  description = NA_character_,
  packages = character(0),
  usesResolution = TRUE
)
```

## Arguments

- fn:

  A function with signature `function(snn, resolution, verbose, ...)`
  returning one cluster id per column of `snn` (vector, factor, or
  single-column data.frame/matrix). See details.

- name:

  Short name for the method, used in messages and in the `algorithm`
  argument of
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md).

- description:

  One-line description of what the method does.

- packages:

  Character vector of packages the method needs. They are checked with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) before the
  method runs, so a missing dependency produces a clear error instead of
  a cryptic one.

- usesResolution:

  Whether the method's output actually depends on `resolution`. Set to
  `FALSE` for algorithms that determine the number of clusters
  themselves (e.g. affinity propagation);
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
  then runs the method once rather than repeating identical work for
  every requested resolution. Defaults to `TRUE`.

## Value

An object of class `gsClusterMethod`.

## Details

`snn` is a sparse symmetric matrix (a Seurat `Graph`) whose rows/columns
are gene sets and whose values are SNN similarities in `[0, 1]`.
Returned ids may be 0- or 1-based;
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
renumbers them to start at 0 and orders the factor levels numerically.
If the returned ids are named, they are reordered to match
`colnames(snn)`.

## See also

[`registerClusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/registerClusterMethod.md),
[`listClusterMethods()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/listClusterMethods.md),
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)

## Examples

``` r
# A backend built on igraph's fast-greedy modularity optimization
fastGreedy <- clusterMethod(
  fn = function(snn, resolution, verbose = FALSE, ...) {
    g <- snnToIgraph(snn)
    igraph::membership(
      igraph::cluster_fast_greedy(g, weights = igraph::E(g)$weight)
    )
  },
  name = "fast_greedy",
  description = "igraph fast-greedy modularity (ignores resolution)",
  packages = "igraph",
  usesResolution = FALSE
)
fastGreedy
#> <gsClusterMethod> fast_greedy
#>   igraph fast-greedy modularity (ignores resolution)
#>   requires: igraph
#>   picks the number of clusters itself (resolution is ignored)

# Try it on the SNN graph of the example data
ids <- fastGreedy$fn(exampleGsList@cluster@graphs$GSEA_snn, resolution = 1)
table(ids)
#> ids
#>  1  2  3  4 
#> 81 19 27 70 

# Use it in the workflow
geneSetsList <- initializeList(exampleGeneSets, exampleEnrichment,
                               verbose = FALSE)
gsListObject <- clusterGeneSets(geneSetsList, algorithm = fastGreedy,
                                verbose = FALSE)
```
