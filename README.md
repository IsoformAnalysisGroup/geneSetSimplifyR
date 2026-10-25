# geneSetSimplifyR

Gene-set enrichment analysis often returns hundreds or thousands of
significant gene sets, far too many to read one by one. `geneSetSimplifyR`
groups them into a small number of clusters of overlapping gene sets and
gives each cluster an informative two-word label.

- **Clustering.** Gene sets are compared by the genes they share, using
  latent semantic indexing and a shared nearest neighbor graph. Several
  community-detection methods are available, and you can plug in your own.
- **Labeling.** Each cluster is named after the word pair that best describes
  its gene-set names (TF-IDF weighted by how central the gene set is).
- **Visualization.** UMAP, bar chart, lollipop, violin, heatmap, dendrogram
  and clustering-tree plots, with support for results from several
  conditions or analyses.

It works on the output of any enrichment tool, with helpers for
[fgsea](https://bioconductor.org/packages/fgsea/) and pairedGSEA results.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("remotes") # needed to install from GitHub
BiocManager::install("IsoformAnalysisGroup/geneSetSimplifyR")
```

`BiocManager` installs the dependencies from both CRAN and Bioconductor.

### Platform notes

- **Linux:** igraph needs the GLPK library. On Debian/Ubuntu, run
  `sudo apt-get install libglpk40` first. For other missing system libraries,
  `pak::pkg_install("IsoformAnalysisGroup/geneSetSimplifyR")` installs the
  package together with the system libraries it needs.
- **macOS:** if installation stops with `ERROR: sub-architecture '...' is not
  installed`, add `INSTALL_opts = "--no-multiarch"` to the
  `BiocManager::install()` call.

## Quick start

```r
library(geneSetSimplifyR)

# A named list of gene sets and, optionally, a table of enrichment results
# with the columns 'pathway', 'enrichment_score' (or 'NES'), and optionally
# 'pval' and 'source'.
gsList <- geneSetSimplifyR(
  geneSetsList = exampleGeneSets,
  geneSetsDF = exampleEnrichment
)

gsList
plotUMAP(gsList)
plotEnrichments(gsList)

# One row per gene set with its cluster label, for supplementary tables
extractTable(gsList, summarize = FALSE)
```

## Documentation

The full documentation, including the vignette and all help pages, is online
at <https://isoformanalysisgroup.github.io/geneSetSimplifyR/>. The vignette
walks through the whole workflow, from fgsea results to renaming clusters:
[Get started](https://isoformanalysisgroup.github.io/geneSetSimplifyR/articles/geneSetSimplifyR.html).

GitHub installs skip the vignette by default. To also read it offline in R,
install with `build_vignettes = TRUE`. This needs a few extra packages and
[pandoc](https://pandoc.org/installing.html), which comes with RStudio:

```r
BiocManager::install(c("knitr", "rmarkdown", "BiocStyle"))
BiocManager::install(
  "IsoformAnalysisGroup/geneSetSimplifyR",
  build_vignettes = TRUE,
  force = TRUE # reinstall if the package is already installed
)
vignette("geneSetSimplifyR", package = "geneSetSimplifyR")
```

Report bugs and questions at
<https://github.com/IsoformAnalysisGroup/geneSetSimplifyR/issues>.
