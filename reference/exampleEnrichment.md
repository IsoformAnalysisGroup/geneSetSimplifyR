# Example Enrichment Results

fgsea results for the TCGA colon adenocarcinoma data analyzed in
Vitting-Seerup & Sandelin (2017), where gene-set enrichment was done
both for differentially expressed genes ("Expression") and
differentially spliced genes ("Splicing"). Only gene sets with an
adjusted p-value below 0.05 are included, a subset of 197 gene sets.

## Usage

``` r
exampleEnrichment
```

## Format

A tibble with 206 rows (one per gene set and source) and 5 columns:

- source:

  The analysis the gene set was found in: "Expression" or "Splicing".

- pathway:

  The gene set name.

- enrichment_score:

  Normalized enrichment score (NES) from fgsea.

- pval:

  p-value from fgsea.

- padj:

  Benjamini-Hochberg adjusted p-value from fgsea.

## Source

Computed with
[`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html)
on MSigDB gene sets; see `inst/scripts/make-exampleEnrichment.R`.
Vitting-Seerup K, Sandelin A. The Landscape of Isoform Switches in Human
Cancers. Mol Cancer Res. 2017;15(9):1206-1220.
