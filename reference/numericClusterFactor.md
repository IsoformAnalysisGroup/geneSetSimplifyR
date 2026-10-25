# Coerce a cluster-id vector/factor to consecutive ids starting at 0

Internal helper. Backends differ in how they number clusters: Seurat is
0-based, igraph is 1-based, a custom backend may return gaps (e.g.
`c(1, 3, 5)`) or non-numeric labels. Seurat may also store factor levels
in lexicographic order ("10" before "2"), which breaks code that pairs
labels to levels by position. This maps the ids to consecutive integers
starting at 0, in numeric order when every id is numeric (character
order otherwise), and returns a factor with numeric level order.

## Usage

``` r
numericClusterFactor(x)
```

## Arguments

- x:

  A factor or vector of cluster ids.

## Value

A factor with levels `0, 1, ..., k - 1`.
