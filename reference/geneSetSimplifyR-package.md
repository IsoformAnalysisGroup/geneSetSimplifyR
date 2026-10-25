# geneSetSimplifyR: Simplify Gene-Set Enrichment Analysis

Simplifies the results of gene-set enrichment analysis (GSEA) by
clustering the enriched gene sets on the genes they share. Gene sets are
embedded with latent semantic indexing (LSI) of their gene membership,
and clustered by community detection on a shared nearest neighbor (SNN)
graph. Each cluster is labeled with the two-word phrase that best
describes its gene-set names, found by term frequency-inverse document
frequency (TF-IDF) and gene-set centrality. Results from several
conditions or analyses can be combined. The package contains several
plotting functions so the user can quickly visualize the clustering
output and draw biologically meaningful insights.

## See also

Useful links:

- <https://github.com/IsoformAnalysisGroup/geneSetSimplifyR>

- Report bugs at
  <https://github.com/IsoformAnalysisGroup/geneSetSimplifyR/issues>

## Author

**Maintainer**: Kristoffer Vitting-Seerup <krivi@dtu.dk>
([ORCID](https://orcid.org/0000-0002-6450-0608))

Authors:

- Ella Hedeboe <ellahedeboe@hotmail.com>
  ([ORCID](https://orcid.org/0000-0002-3978-4661))

- Chiao-Yu Hsieh <chihs@dtu.dk>
  ([ORCID](https://orcid.org/0009-0007-0463-5977))
