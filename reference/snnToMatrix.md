# Convert an SNN graph to a dense similarity matrix

Helper for writing clustering methods that operate on a full similarity
matrix rather than a graph (affinity propagation, binary cut, k-means on
similarity profiles). The SNN graph is already a similarity matrix with
values in `[0, 1]`; this simply densifies it and keeps the gene set
names as dimnames.

## Usage

``` r
snnToMatrix(snn, warnAbove = 5000)
```

## Arguments

- snn:

  The SNN graph passed to a cluster method.

- warnAbove:

  Number of gene sets above which a warning about the memory cost of
  densifying is emitted. Set to `Inf` to silence. The dense matrix needs
  roughly `8 * N^2` bytes, so 10,000 gene sets is about 800 MB.

## Value

A dense symmetric numeric matrix with gene sets as rows/columns.

## See also

[`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md),
[`snnToIgraph()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/snnToIgraph.md)

## Examples

``` r
simMat <- snnToMatrix(exampleGsList@cluster@graphs$GSEA_snn)
simMat[1:3, 1:3]
#>                                           GRADE_COLON_AND_RECTAL_CANCER_UP
#> GRADE_COLON_AND_RECTAL_CANCER_UP                                         1
#> BENPORATH_NANOG_TARGETS                                                  0
#> GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP                                0
#>                                           BENPORATH_NANOG_TARGETS
#> GRADE_COLON_AND_RECTAL_CANCER_UP                                0
#> BENPORATH_NANOG_TARGETS                                         1
#> GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP                       0
#>                                           GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP
#> GRADE_COLON_AND_RECTAL_CANCER_UP                                                  0
#> BENPORATH_NANOG_TARGETS                                                           0
#> GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP                                         1
```
