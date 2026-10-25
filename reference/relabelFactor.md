# Relabel the levels of a factor

Internal helper replacing
[`plyr::mapvalues()`](https://rdrr.io/pkg/plyr/man/mapvalues.html) for
factors: every level found in `from` is renamed to the matching element
of `to`; other levels are kept.

## Usage

``` r
relabelFactor(x, from, to)
```

## Arguments

- x:

  A factor.

- from, to:

  Character vectors of equal length.

## Value

The relabeled factor.
