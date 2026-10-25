# Extract summary table of clustering

Extract summary table of clustering

## Usage

``` r
extractTable(
  geneSetsList,
  resolution = NULL,
  summarize = TRUE,
  removeClusterId = TRUE
)
```

## Arguments

- geneSetsList:

  'gsList' object with clusters and labels.

- resolution:

  Numeric. Resolution to extract the table for. Default `NULL` uses the
  active resolution: the optimal resolution found by
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md),
  or the one last set with
  [`updateResolution()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateResolution.md)
  or
  [`updateClusterLabels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateClusterLabels.md).

- summarize:

  Logical indicating whether to summarize the table by cluster (and
  source). Default is TRUE.

- removeClusterId:

  Logical, if TRUE the cluster number will be removed from each cluster
  label. Default: TRUE.

## Value

A data.frame. With `summarize = TRUE`, one row per cluster and source
with the mean and standard deviation of the enrichment score (ignoring
missing values) and the number of gene sets. With `summarize = FALSE`,
one row per gene set and source with all columns of the
'enrichmentScore' slot plus its `cluster`; if a cluster label was
changed with
[`updateClusterLabels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateClusterLabels.md),
the original label is given in `initial_cluster_label`.

## Examples

``` r
extractTable(exampleGsList)
#>       source               cluster enrichment mean enrichment sd genesets
#> 1 Expression             Stem Cell            1.27          0.07       81
#> 2   Splicing             Stem Cell            1.14            NA        1
#> 3 Expression        Serum Response            1.21          0.08       31
#> 4   Splicing        Serum Response            1.19          0.05       44
#> 5 Expression       Peptidyl Lysine            1.25            NA        1
#> 6   Splicing       Peptidyl Lysine            1.23          0.06       25
#> 7 Expression Carbonate Dehydratase            1.21          0.13       15
#> 8   Splicing Carbonate Dehydratase            1.23          0.17        8
head(extractTable(exampleGsList, resolution = 0.6, summarize = FALSE))
#>       source                                   pathway enrichment_score
#> 1 Expression GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP         1.294445
#> 2 Expression                     FISCHER_DREAM_TARGETS         1.223318
#> 3 Expression         MEBARKI_HCC_PROGENITOR_FZD8CRD_UP         1.282581
#> 4 Expression            VECCHI_GASTRIC_CANCER_EARLY_UP         1.283857
#> 5 Expression    SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP         1.420683
#> 6 Expression          KOBAYASHI_EGFR_SIGNALING_24HR_DN         1.340100
#>           pval         padj nlog10pval   cluster
#> 1 3.637638e-15 2.210229e-11  14.439181 Stem Cell
#> 2 3.314153e-14 1.006840e-10  13.479627 Stem Cell
#> 3 3.863003e-13 9.388643e-10  12.413075 Stem Cell
#> 4 2.789000e-11 4.236490e-08  10.554552 Stem Cell
#> 5 1.832008e-10 2.023870e-07   9.737073 Stem Cell
#> 6 3.091404e-10 2.889749e-07   9.509844 Stem Cell
```
