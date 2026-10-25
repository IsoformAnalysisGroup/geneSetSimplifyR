# Resolve which resolution a function should use

Internal helper shared by the plotting and table functions. An explicit
`resolution` always wins. With `resolution = NULL` the active resolution
is used: the one last set by
[`updateClusterLabels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateClusterLabels.md)
or
[`updateResolution()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateResolution.md),
or the optimal resolution found by
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
if neither was used.

## Usage

``` r
resolveResolutionColumn(geneSetsList, resolution = NULL)
```

## Arguments

- geneSetsList:

  A clustered 'gsList' object.

- resolution:

  `NULL` or a single resolution that was clustered.

## Value

The name of the resolution column, e.g. `"GSEA_snn_res.1.1"`.
