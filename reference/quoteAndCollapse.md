# Quote and comma-collapse a character vector for use in messages

Internal helper.

## Usage

``` r
quoteAndCollapse(x)
```

## Arguments

- x:

  A character vector.

## Value

A single string, e.g. `"'pval', 'padj'"`. Returns `"none"` for an empty
vector so error messages stay readable.
