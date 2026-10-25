# Update default resolution of gsList

Sets the resolution used by all plotting and table functions when their
`resolution` argument is `NULL`. It replaces the optimal resolution
found by
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
as the default.

## Usage

``` r
updateResolution(geneSetsList, resolution)
```

## Arguments

- geneSetsList:

  'gsList' object that has been clustered and labeled.

- resolution:

  Numeric. The clustering resolution to use as the new default. It must
  have been used as a resolution in
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md).

## Value

'gsList' object updated with the new default resolution.

## Examples

``` r
gsListObject <- updateResolution(exampleGsList, 0.6)
extractResolution(gsListObject)
#> [1] 0.6
```
