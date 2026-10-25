# Rename a cluster

Replaces the label of one cluster with a name of your choice. The
original label is kept and reported as `initial_cluster_label` by
`extractTable(summarize = FALSE)`, also after several renames.

## Usage

``` r
updateClusterLabels(
  geneSetsList,
  resolution = NULL,
  clusterNumber,
  newClusterName
)
```

## Arguments

- geneSetsList:

  'gsList' object that has been clustered and labeled.

- resolution:

  Numeric. The clustering resolution of the cluster to rename. It must
  have been used as a resolution in
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md).
  When given, this resolution also becomes the active resolution:
  plotting and table functions called with `resolution = NULL`
  afterwards use it. Default `NULL` renames a cluster of the active
  resolution.

- clusterNumber:

  The number of the cluster to rename, as shown in front of its label
  (e.g. `0` for "0: Immune Response").

- newClusterName:

  The new name of the cluster.

## Value

'gsList' object with the cluster label updated.

## Examples

``` r
gsListObject <- updateClusterLabels(
  geneSetsList = exampleGsList,
  clusterNumber = 0,
  newClusterName = 'cancer response'
)
head(extractTable(gsListObject, summarize = FALSE))
#>       source                                   pathway enrichment_score
#> 1 Expression GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP         1.294445
#> 2 Expression                     FISCHER_DREAM_TARGETS         1.223318
#> 3 Expression         MEBARKI_HCC_PROGENITOR_FZD8CRD_UP         1.282581
#> 4 Expression            VECCHI_GASTRIC_CANCER_EARLY_UP         1.283857
#> 5 Expression    SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP         1.420683
#> 6 Expression          KOBAYASHI_EGFR_SIGNALING_24HR_DN         1.340100
#>           pval         padj nlog10pval         cluster initial_cluster_label
#> 1 3.637638e-15 2.210229e-11  14.439181 cancer response             Stem Cell
#> 2 3.314153e-14 1.006840e-10  13.479627 cancer response             Stem Cell
#> 3 3.863003e-13 9.388643e-10  12.413075 cancer response             Stem Cell
#> 4 2.789000e-11 4.236490e-08  10.554552 cancer response             Stem Cell
#> 5 1.832008e-10 2.023870e-07   9.737073 cancer response             Stem Cell
#> 6 3.091404e-10 2.889749e-07   9.509844 cancer response             Stem Cell
```
