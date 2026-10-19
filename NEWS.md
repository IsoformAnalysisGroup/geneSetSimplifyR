# geneSetSimplifyR 0.99.0

First version submitted to Bioconductor.

## New features

* Pluggable clustering methods: `clusterGeneSets(algorithm =)` accepts a
  registered name ('louvain', 'louvain.multilevel', 'slm', 'leiden',
  'igraph_louvain', 'igraph_leiden', 'kmeans', 'apcluster', 'binary_cut'), a
  method made with `clusterMethod()`, or a function. See
  `registerClusterMethod()`, `listClusterMethods()`, `getClusterMethod()`,
  `snnToIgraph()` and `snnToMatrix()`.
* `clusterGeneSets(K =)` sets the number of nearest neighbors; the default is
  `round(sqrt(N))` for `N` gene sets.
* Parsers for enrichment results: `parseFgsea()` (GSEA and ORA results) and
  `parsePairedGSEA()`. `getGeneSets()` fetches MSigDB gene sets.
* `initializeList()` accepts a `pval` (or `pvalue`) column and adds
  `nlog10pval`. Other numeric columns are kept.
* `initializeList()` and `geneSetSimplifyR()` accept the path to a GMT file
  as `geneSetsList`.
* `plotEnrichments()`, `plotHeatmap()` and `plotViolin()` gain `plotFrom` (the
  column to plot; by default -log10 p-values when available, else enrichment
  scores). `plotEnrichments()` and `plotHeatmap()` gain `summaryStat`.
* `plotHeatmap(scaleCols =)` with default `"auto"`: clusters are z-scored
  across sources only when there are three or more sources.
* `labelClusters(removeFirstWord =)` drops the first word (e.g. an author
  name or database prefix) of gene-set names before labeling.
* `geneSetSimplifyR()` passes `K`, `stopWords` and `nLabels` on.
* `normalizeClusterLevels()` updates objects clustered by earlier versions.
* The `gsList` class has a `metadata` slot recording parameters and results.

## Changes in behavior

* With `resolution = NULL`, plotting and table functions use the *active*
  resolution: the optimal one from `clusterGeneSets()`, or the one last set
  with `updateResolution()` or `updateClusterLabels(resolution =)`.
* The LSI is computed internally (same results as before) instead of with
  Signac, which removes Signac and its genomics dependencies.
* TF-IDF method 1 is used (chosen by grid search).
* Clusters always get distinct labels; if needed the cluster id is appended.
* `plotDendrogram()` returns a ggplot instead of drawing a base plot.
* `plotHeatmap()` returns the pheatmap object invisibly and shows clusters
  absent from a source in gray instead of as 0.
* `plotProportions()` counts a gene set found in several sources once, as
  "Both" or "Multiple", and requires at least two sources.
* Re-clustering a labeled object removes its old labels.
* `clusterGeneSets()` no longer changes the random number stream of the
  session.

## Bug fixes

* `plotUMAP(interactive = TRUE)` works.
* `plotHeatmap()` works with a single source and with constant values.
* `updateClusterLabels()` keeps the original label when called repeatedly,
  validates its input, and `extractTable(summarize = FALSE)` reports it as
  `initial_cluster_label`.
* `plotViolin(removeClusterId = TRUE)` keeps the clusters ordered by median.
* Inputs with 50 or fewer gene sets can be clustered (at least 10 needed).
* A resolution where the silhouette is undefined no longer stops
  `clusterGeneSets()`; if all resolutions fail, the error lists why.
* The LSI falls back to another SVD solver when irlba does not converge.
* Stop words with special characters, and an empty stop-word list, work.
* Long gene-set names keep their "_UP"/"_DN" suffix when shortened.
* `plotClustree()` no longer attaches ggraph to the search path.
* Clear errors for unknown resolutions and cluster numbers.

# geneSetSimplifyR 1.0.0 (2026-05)

* Internal release, not published.
* Updated the code to adapt to the latest version of Seurat.
* Renamed `highligthClusterNo` to `highlightClusterNo` (and its
  `highligthPoint*` / `highligthLabel*` siblings) in all plotting functions.
* Renamed `plotDendogram()` to `plotDendrogram()`.
* `weighted_score` in label generation weighs TF-IDF and centrality scores
  equally.
* Added `LazyDataCompression: xz`, reducing package size.

# geneSetSimplifyR 0.1.3 (2023-12)

* Updates to fix errors introduced by a Matrix update.

# geneSetSimplifyR 0.1.2 (2023-06)

* Updated description.
* `clusterGeneSets()` has new (slightly changed) default parameters based on
  a new grid search.
* `clusterGeneSets()` detects the best clustering resolution (of the ones
  tested) and sets it as the default used throughout the package.
* Plot functions handle user-supplied resolutions (`plotUMAP`,
  `plotProportions`, `plotViolin`, `plotClustree`) and different numbers of
  conditions (`plotUMAP`, `plotViolin`).
* New `plotClusterLabels()` and `updateClusterLabels()` for in-depth analysis
  and fine-tuning of cluster labels.
* `plotUMAP`, `plotProportions` and `plotViolin` can highlight clusters.
* All functions use `resolution` as argument instead of `defaultResolution`.
* Added a show method for gsList that summarizes the object.

# geneSetSimplifyR 0.1.0 (2022-10-10)

* Package established.
