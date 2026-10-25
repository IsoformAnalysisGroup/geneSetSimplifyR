# The gsList class

A `gsList` holds the input and results of the geneSetSimplifyR workflow.
It is created by
[`initializeList()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/initializeList.md),
filled by
[`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
and
[`labelClusters()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/labelClusters.md),
or all at once by
[`geneSetSimplifyR()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/geneSetSimplifyR.md).
Use
[`extractTable()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/extractTable.md)
to get the results as a table.

## Usage

``` r
# S4 method for class 'gsList'
show(object)
```

## Arguments

- object:

  A 'gsList' object.

## Value

`gsList` objects are created by
[`initializeList()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/initializeList.md),
not by calling the class generator directly. `show()` prints a summary
of the object and returns it invisibly.

## Functions

- `show(gsList)`: Print a short summary of the object.

## Slots

- `genesetsDataframe`:

  A data.frame with genes as rows and gene sets as columns; 1 if the
  gene is in the gene set, else 0.

- `enrichmentScore`:

  A data.frame with one row per gene set and source: `pathway`,
  `source`, `enrichment_score`, `pval`, `nlog10pval` and any other
  column supplied to
  [`initializeList()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/initializeList.md).

- `cluster`:

  A Seurat object with the LSI, UMAP, SNN graph and the cluster of every
  gene set for each resolution (`GSEA_snn_res.<resolution>` metadata
  columns).

- `labels`:

  A list with, for each resolution, the cluster labels
  (`labelsDefault`), the TF-IDF scores of the label candidates
  (`labelsTFIDF`) and the most central gene sets (`labelsCentrality`).

- `metadata`:

  A list with parameters and results of the workflow, e.g. the number of
  nearest neighbors, the chosen resolution and silhouette scores.

## Examples

``` r
exampleGsList
#> A gsList object with 197 gene sets from 2 source(s) clustered into 4 labeled clusters.
#> Metadata:
#>   algorithm: louvain
#>   k (KNN): 14
#>   LSI elbow cutoff: 11
#>   selected resolution: 0.1
#>   n pathways: 197
names(exampleGsList@metadata)
#>  [1] "n_genes"             "n_pathways"          "n_sources"          
#>  [4] "package_version"     "initialized_at"      "k"                  
#>  [7] "lsi_elbow_cutoff"    "cluster_algorithm"   "resolutions_tested" 
#> [10] "selected_resolution" "active_resolution"   "silhouette_scores"  
#> [13] "clustered_at"        "n_labels"            "remove_first_word"  
#> [16] "stop_words"          "labeled_at"         
```
