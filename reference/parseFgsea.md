# Parse fgsea results into geneSetSimplifyR input

Converts [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html)
output into the two inputs
[`geneSetSimplifyR()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/geneSetSimplifyR.md)
needs: `geneSetsList` and `geneSetsDF` (with `source`, `pathway`,
`enrichment_score`). To combine several fgsea runs into one multi-source
analysis, call this function once per run with a distinct `source` label
and [`rbind()`](https://rdrr.io/r/base/cbind.html) the resulting
`geneSetsDF`s.

## Usage

``` r
parseFgsea(fgseaResult, geneSets, cutoff = 0.05, source = "fgsea")
```

## Arguments

- fgseaResult:

  Output of [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html)
  or [`fgsea::fora()`](https://rdrr.io/pkg/fgsea/man/fora.html). Must
  contain `pathway` and `padj` columns, and either `NES` (GSEA) or
  `foldEnrichment` (over-representation analysis, used as
  `log2(foldEnrichment)`).

- geneSets:

  Named list of gene-set memberships (the same one passed to
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html)).

- cutoff:

  padj threshold for filtering pathways. Default 0.05.

- source:

  Label written to the `source` column of `geneSetsDF`. Use a distinct
  value per run when combining multiple analyses. Default `"fgsea"`.

## Value

A list with:

- `geneSetsList`:

  `geneSets` filtered to pathways passing the cutoff.

- `geneSetsDF`:

  Long-format data frame: `source`, `pathway`, `enrichment_score`
  (fgsea's NES), `pval` (when present) and `padj`. One row per pathway.

## Examples

``` r
# A small table in the format of fgsea::fgsea() output
pathways <- names(exampleGeneSets)[1:60]
fgseaRes <- data.frame(
  pathway = pathways,
  pval = rep(c(0.001, 0.01, 0.3), 20),
  padj = rep(c(0.004, 0.03, 0.5), 20),
  NES = seq(-2.5, 2.5, length.out = 60)
)
parsed <- parseFgsea(
  fgseaResult = fgseaRes,
  geneSets    = exampleGeneSets,
  cutoff      = 0.05,
  source      = "treatmentA"
)
head(parsed$geneSetsDF)
#>       source                                   pathway enrichment_score  pval
#> 1 treatmentA GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP        -2.500000 0.001
#> 2 treatmentA                     FISCHER_DREAM_TARGETS        -2.415254 0.010
#> 3 treatmentA                       LEE_BMP2_TARGETS_DN        -2.245763 0.001
#> 4 treatmentA            VECCHI_GASTRIC_CANCER_EARLY_UP        -2.161017 0.010
#> 5 treatmentA    SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP        -1.991525 0.001
#> 6 treatmentA           CAIRO_HEPATOBLASTOMA_CLASSES_UP        -1.906780 0.010
#>    padj
#> 1 0.004
#> 2 0.030
#> 3 0.004
#> 4 0.030
#> 5 0.004
#> 6 0.030
gss <- geneSetSimplifyR(
  geneSetsList = parsed$geneSetsList,
  geneSetsDF   = parsed$geneSetsDF,
  resolution   = 1.1,
  verbose      = FALSE
)
```
