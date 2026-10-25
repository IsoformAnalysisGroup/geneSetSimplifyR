# Run a clustering method on an SNN graph

Internal. Checks the method's dependencies, calls it for a single
resolution, and validates/orders the returned cluster ids.

## Usage

``` r
runClusterMethod(
  method,
  snn,
  resolution,
  verbose = FALSE,
  clusterArgs = list()
)
```

## Arguments

- method:

  A `gsClusterMethod`.

- snn:

  The SNN graph (sparse symmetric matrix).

- resolution:

  A single numeric resolution.

- verbose:

  Passed to the method.

- clusterArgs:

  Named list of extra arguments forwarded to the method.

## Value

An unnamed vector of cluster ids, ordered as `colnames(snn)`.
