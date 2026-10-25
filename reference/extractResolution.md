# Extract the resolution currently used by default

Returns the resolution that the plotting and table functions use when
their `resolution` argument is `NULL`. This is the optimal resolution
found by
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md),
unless it was changed with
[`updateResolution()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateResolution.md)
or `updateClusterLabels(resolution = ...)`.

## Usage

``` r
extractResolution(geneSetsList)
```

## Arguments

- geneSetsList:

  'gsList' object that has been clustered with the 'clusterGeneSets()'
  function.

## Value

The default resolution, as a number.

## Examples

``` r
extractResolution(exampleGsList)
#> [1] 0.1
```
