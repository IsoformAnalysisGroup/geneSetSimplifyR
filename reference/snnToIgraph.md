# Convert an SNN graph to an igraph object

Helper for writing igraph-based clustering methods. Turns the sparse
symmetric SNN matrix handed to a cluster method into a weighted
undirected `igraph` graph whose vertex names are the gene sets.

## Usage

``` r
snnToIgraph(snn)
```

## Arguments

- snn:

  The SNN graph passed to a cluster method (a sparse symmetric matrix /
  Seurat `Graph`).

## Value

A weighted undirected `igraph` object.

## See also

[`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md)

## Examples

``` r
g <- snnToIgraph(exampleGsList@cluster@graphs$GSEA_snn)
g
#> IGRAPH bd6aa5a UNW- 197 3750 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from bd6aa5a (vertex names):
#> [1] GRADE_COLON_AND_RECTAL_CANCER_UP--LEE_BMP2_TARGETS_DN                                 
#> [2] GRADE_COLON_AND_RECTAL_CANCER_UP--RODRIGUES_THYROID_CARCINOMA_POORLY_DIFFERENTIATED_UP
#> [3] GRADE_COLON_AND_RECTAL_CANCER_UP--CAIRO_HEPATOBLASTOMA_CLASSES_UP                     
#> [4] GRADE_COLON_AND_RECTAL_CANCER_UP--BERENJENO_TRANSFORMED_BY_RHOA_UP                    
#> [5] GRADE_COLON_AND_RECTAL_CANCER_UP--MANALO_HYPOXIA_DN                                   
#> [6] GRADE_COLON_AND_RECTAL_CANCER_UP--WANG_RESPONSE_TO_GSK3_INHIBITOR_SB216763_DN         
#> + ... omitted several edges
```
