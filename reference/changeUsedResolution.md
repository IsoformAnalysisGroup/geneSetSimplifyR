# Set the cluster identities used for plotting

Internal helper. Sets `active.ident` of the Seurat object to the
clusters of the requested resolution, labeled as `"<id>: <label>"` (or
`"<label>"` with `removeClusterId = TRUE`). Only the returned copy is
changed, so plotting a different resolution never alters the user's
object.

## Usage

``` r
changeUsedResolution(geneSetsList, resolution = NULL, removeClusterId = FALSE)
```

## Arguments

- geneSetsList:

  A labeled 'gsList' object.

- resolution:

  `NULL` (the active resolution) or a clustered resolution.

- removeClusterId:

  Logical. Drop the `"<id>: "` prefix from the labels.

## Value

The 'gsList' object with updated identities.
