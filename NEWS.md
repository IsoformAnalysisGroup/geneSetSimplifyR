# geneSetSimplifyR 0.99.0

* First version submitted to Bioconductor.
* `geneSetSimplifyR()` clusters and labels enriched gene sets in one call. The
  steps are also available on their own: `initializeList()`,
  `clusterGeneSets()` and `labelClusters()`. Clustering needs more than 50
  gene sets.
* Gene sets are clustered by community detection on a shared nearest neighbor
  (SNN) graph built from their gene membership. `clusterGeneSets(algorithm =)`
  accepts a registered method ('louvain', 'louvain.multilevel', 'slm',
  'leiden', 'igraph_louvain', 'igraph_leiden', 'kmeans', 'apcluster',
  'binary_cut'), a method made with `clusterMethod()`, or a function. See
  `registerClusterMethod()`, `listClusterMethods()`, `getClusterMethod()`,
  `getSNN()`, `snnToIgraph()` and `snnToMatrix()`.
* Clusters are labeled with the two-word phrase that best describes their
  gene-set names. Words in `stopWords` are ignored, and
  `labelClusters(removeFirstWord =)` drops author names or database prefixes.
* Several resolutions are clustered at once and the one with the highest
  silhouette score is used by default. Change it with `updateResolution()`
  and rename clusters with `updateClusterLabels()`.
* Plots: `plotBarchart()`, `plotUMAP()`, `plotEnrichments()`, `plotHeatmap()`,
  `plotViolin()`, `plotProportions()`, `plotDendrogram()`, `plotClustree()`
  and `plotClusterLabels()`. `extractTable()` returns the results as a table.
* `parseFgsea()` (GSEA and ORA results) and `parsePairedGSEA()` turn
  enrichment results into input for `geneSetSimplifyR()`. `getGeneSets()`
  fetches MSigDB gene sets, and gene sets can also be given as the path to a
  GMT file.
* Example data: `exampleGeneSets`, `exampleEnrichment` and `exampleGsList`.
