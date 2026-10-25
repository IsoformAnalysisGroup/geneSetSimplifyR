# Determine the elbow cutoff to be used to run UMAP

Determine the elbow cutoff to be used to run UMAP

## Usage

``` r
kvsElbow(variance, method = c("distance", "zscore"), zScoreCutoff = 2.5)
```

## Arguments

- variance:

  Standard deviations of the LSI components, stored in a Seurat object
  as `seuratObj@reductions$lsi@stdev`. Must be sorted in decreasing
  order.

- method:

  Which elbow-detection method to use. One of:

  `"distance"` (default)

  :   Returns the point of the (smoothed) scree curve that lies farthest
      from the chord connecting its endpoints. Light spline smoothing
      keeps it from being thrown off by jitter in the noisy tail, and it
      needs no threshold.

  `"zscore"`

  :   Flags the stepwise drops (1st derivative) that are statistical
      outliers via a robust z-score and returns the last such drop.
      Sensitive to the `zScoreCutoff` threshold.

- zScoreCutoff:

  Integer indicating the cut off value for the Z-score. Only used when
  `method = "zscore"`.

## Value

A list with element `elbow_cutoff`: the index of the last LSI component
to use (at least 2).

## Details

Short scree curves (fewer than 50 components, i.e. inputs with 50 or
fewer gene sets) are supported. With fewer than four components after
dropping the first one, all components are used.
