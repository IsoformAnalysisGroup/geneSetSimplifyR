# Build display labels aligned to a set of cluster ids

Internal helper. Takes the `labelsDefault` table for one resolution,
keeps the top-scoring label per cluster, and returns one label per
element of `clusterIds`, in that order.

## Usage

``` r
alignedClusterLabels(
  geneSetsList,
  resolutionColumn,
  clusterIds,
  removeClusterId = FALSE
)
```

## Arguments

- geneSetsList:

  A labeled 'gsList' object.

- resolutionColumn:

  Name of the resolution, e.g. `"GSEA_snn_res.1"`.

- clusterIds:

  Character vector of cluster ids to label, in the desired output order
  (typically `levels(active.ident)`).

- removeClusterId:

  Logical. Drop the `"<id>: "` prefix from the labels.

## Value

Character vector of labels, one per element of `clusterIds`.

## Details

Three things this guards against:

- **Order.** The label table is in numeric cluster order; a factor's
  levels need not be, so labels are matched by cluster id rather than by
  position.

- **Coverage.** A cluster can be missing from the label table (e.g. a
  cluster whose gene-set names yield no usable bigram). Those fall back
  to the bare cluster id, so the returned vector always has the same
  length as `clusterIds`.

- **Uniqueness.** When `removeClusterId = TRUE`, two clusters with the
  same label (e.g. after a manual rename) would otherwise be merged into
  a single factor level. Those clusters keep their id prefix instead.
