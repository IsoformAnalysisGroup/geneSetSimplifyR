# Cluster gene sets

Cluster gene sets

## Usage

``` r
clusterGeneSets(
  geneSetsList,
  resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 4.1, 5.1),
  algorithm = "louvain",
  clusterArgs = list(),
  K = NA,
  elbow_method = "distance",
  verbose = TRUE
)
```

## Arguments

- geneSetsList:

  'gsList' object created with the 'initializeList()' function.

- resolution:

  Vector with the desired clustering resolutions to use. Higher number
  usually results in more clusters.

- algorithm:

  Clustering method applied to the SNN graph. Either the name of a
  registered method (see
  [`listClusterMethods()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/listClusterMethods.md);
  e.g. 'louvain', 'leiden', 'slm', 'igraph_louvain', 'igraph_leiden',
  'kmeans', 'apcluster', 'binary_cut'), a `gsClusterMethod` built with
  [`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md),
  or a function following the same contract. Only this step is swappable
  – the SNN graph itself is always built the same way.

  Methods that determine the number of clusters themselves ('kmeans',
  'apcluster', 'binary_cut') ignore `resolution` and are run once; the
  result is stored under the first requested resolution.

- clusterArgs:

  Named list of extra arguments passed on to the clustering method, for
  example `list(n.start = 20)` for the Seurat backends,
  `list(maxK = 50)` for 'kmeans', `list(p = -1)` for 'apcluster', or
  `list(cutoff = 0.9)` for 'binary_cut'.

- K:

  Integer or NA. Number of nearest neighbors used when building the KNN
  graph. If `NA` (default) it is set to `round(sqrt(N))`, where `N` is
  the number of gene sets. Must be smaller than `N`.

- elbow_method:

  Method for choosing the number of LSI dimensions from the scree curve:
  'distance' (default) or 'zscore'. See the Details.

- verbose:

  a logic indicating whether to print progress statements. Defaults to
  TRUE.

## Value

'gsList' object with clustering stored. The best resolution (highest
silhouette score) is annotated and used throughout the rest of the
workflow (unless specified in subsequent functions).

## Details

At least 10 gene sets are needed. The LSI uses up to 50 dimensions
(fewer when there are 50 gene sets or less); the first one, which mostly
reflects gene-set size, is always left out.

Clustering does not change the random number stream of the R session:
the UMAP and the seeded clustering backends restore it afterwards.
Re-clustering an object removes its existing labels, so run
[`labelClusters()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/labelClusters.md)
again afterwards.

## See also

[`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md),
[`registerClusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/registerClusterMethod.md),
[`listClusterMethods()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/listClusterMethods.md)

## Examples

``` r
geneSetsList <- initializeList(
  geneSetsList = exampleGeneSets,
  geneSetsDF = exampleEnrichment,
  verbose = FALSE
)
gsListObject <- clusterGeneSets(
  geneSetsList = geneSetsList,
  resolution = c(0.6, 1.1),
  verbose = FALSE
)

# Swap in a different community-detection implementation
gsListObject <- clusterGeneSets(
  geneSetsList,
  algorithm = "igraph_leiden",
  resolution = 1,
  verbose = FALSE
)
```
