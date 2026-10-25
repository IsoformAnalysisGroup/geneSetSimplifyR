# Line segments of a dendrogram

Line segments of a dendrogram

## Usage

``` r
dendrogramSegments(hc, leafHeight = 0)
```

## Arguments

- hc:

  An `hclust` object.

- leafHeight:

  Height at which the leaves are drawn.

## Value

A data.frame with columns `x`, `y`, `xend` and `yend`; leaves are at
`x = 1, 2, ...` in the order of `hc$order` and at `leafHeight`.
