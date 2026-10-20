## Adds the fgsea p-values to data/exampleEnrichment.rda.
##
## 'exampleEnrichment' is a subset of the gene-set enrichment results for the
## TCGA colon adenocarcinoma (COAD) data analyzed in Vitting-Seerup & Sandelin
## (2017, Molecular Cancer Research 15:1206). Enrichment was computed with
## fgsea::fgseaMultilevel(scoreType = "pos", minSize = 10, maxSize = 1000) on
## the MSigDB gene sets, separately for
##   - differential gene expression ("Expression"), genes ranked by
##     -log10(p-value) * sign(log fold change), and
##   - differential gene usage/splicing ("Splicing"), genes ranked by
##     -log10(p-value).
## The example keeps the results with padj < 0.05 for the gene sets in
## 'exampleGeneSets' (206 gene set/source pairs, 197 gene sets).
##
## Input (not shipped with the package; created by the manuscript analysis):
##   20_1_tcga_paired_fgsea.Rdata with the fgsea results 'tcgaDgeFsc' and
##   'tcgaDguFsc', lists by cancer type.
##
## The rows and their order are kept as in the earlier versions of the
## example, so the clustering of the example is unchanged; the enrichment
## scores are fgsea's NES and are checked against the fgsea results.
##
## Run from the root of the package source:
##   Rscript inst/scripts/make-exampleEnrichment.R path/to/20_1_tcga_paired_fgsea.Rdata
## and then inst/scripts/make-exampleGsList.R.

fgseaFile <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(fgseaFile), file.exists(fgseaFile))

load(fgseaFile) # tcgaDgeFsc, tcgaDguFsc
load(file.path("data", "exampleEnrichment.rda"))
load(file.path("data", "exampleGeneSets.rda"))

fgseaCoad <- rbind(
  data.frame(source = "Expression", as.data.frame(tcgaDgeFsc$COAD)[, c("gene_set", "pval", "padj", "NES")]),
  data.frame(source = "Splicing", as.data.frame(tcgaDguFsc$COAD)[, c("gene_set", "pval", "padj", "NES")])
)

## The example is exactly the significant results for its gene sets
significant <- fgseaCoad[fgseaCoad$padj < 0.05 & !is.na(fgseaCoad$padj) &
                           fgseaCoad$gene_set %in% names(exampleGeneSets), ]
stopifnot(setequal(
  paste(significant$source, significant$gene_set),
  paste(exampleEnrichment$source, exampleEnrichment$pathway)
))

idx <- match(
  paste(exampleEnrichment$source, exampleEnrichment$pathway),
  paste(fgseaCoad$source, fgseaCoad$gene_set)
)
stopifnot(!anyNA(idx))
stopifnot(isTRUE(all.equal(exampleEnrichment$enrichment_score, fgseaCoad$NES[idx])))

exampleEnrichment <- tibble::tibble(
  source = exampleEnrichment$source,
  pathway = exampleEnrichment$pathway,
  enrichment_score = exampleEnrichment$enrichment_score,
  pval = fgseaCoad$pval[idx],
  padj = fgseaCoad$padj[idx]
)

save(
  exampleEnrichment,
  file = file.path("data", "exampleEnrichment.rda"),
  compress = "xz"
)
