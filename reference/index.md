# Package index

## Workflow

Cluster and label gene sets.

- [`geneSetSimplifyR()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/geneSetSimplifyR.md)
  : Run the geneSetSimplifyR workflow
- [`initializeList()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/initializeList.md)
  : Initialize gene sets list
- [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)
  : Cluster gene sets
- [`labelClusters()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/labelClusters.md)
  : Label clusters in clustered gene sets list object

## Visualization

- [`plotBarchart()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotBarchart.md)
  : Plot clusters as a bar chart
- [`plotClusterLabels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotClusterLabels.md)
  : Plot label information for a particular cluster
- [`plotClustree()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotClustree.md)
  : Plot cluster tree
- [`plotDendrogram()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotDendrogram.md)
  : Plot dendrogram of cluster averages
- [`plotEnrichments()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotEnrichments.md)
  : Plot enrichment as lollipop plot
- [`plotHeatmap()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotHeatmap.md)
  : Plot enrichment as heatmap
- [`plotProportions()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotProportions.md)
  : Plot proportions of gene sets in each cluster for each source
- [`plotUMAP()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotUMAP.md)
  : Scatter plot of UMAP dimensionality reduction
- [`plotViolin()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/plotViolin.md)
  : Violin plot of clustering

## Results and resolution

- [`extractTable()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/extractTable.md)
  : Extract summary table of clustering
- [`extractResolution()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/extractResolution.md)
  : Extract the resolution currently used by default
- [`updateResolution()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateResolution.md)
  : Update default resolution of gsList
- [`updateClusterLabels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/updateClusterLabels.md)
  : Rename a cluster
- [`normalizeClusterLevels()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/normalizeClusterLevels.md)
  : Normalize cluster factor level order in a clustered gsList

## Input helpers

Convert enrichment results and fetch gene sets.

- [`parseFgsea()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/parseFgsea.md)
  : Parse fgsea results into geneSetSimplifyR input
- [`parsePairedGSEA()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/parsePairedGSEA.md)
  : Parse pairedGSEA results into geneSetSimplifyR input
- [`getGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/getGeneSets.md)
  : Fetch MSigDB gene sets by pathway name

## Clustering methods

Use or add community-detection methods for the SNN graph.

- [`listClusterMethods()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/listClusterMethods.md)
  : List the available clustering methods

- [`getClusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/getClusterMethod.md)
  : Resolve a clustering method specification

- [`clusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterMethod.md)
  :

  Define a clustering method for
  [`clusterGeneSets()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/clusterGeneSets.md)

- [`registerClusterMethod()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/registerClusterMethod.md)
  : Register a clustering method for use by name

- [`snnToIgraph()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/snnToIgraph.md)
  : Convert an SNN graph to an igraph object

- [`snnToMatrix()`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/snnToMatrix.md)
  : Convert an SNN graph to a dense similarity matrix

## Classes and data

- [`show(`*`<gsList>`*`)`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/gsList-class.md)
  : The gsList class
- [`exampleGeneSets`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/exampleGeneSets.md)
  : Example Gene Sets
- [`exampleEnrichment`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/exampleEnrichment.md)
  : Example Enrichment Results
- [`exampleGsList`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/exampleGsList.md)
  : Example clustered and labeled gene sets
- [`stopWords`](https://isoformanalysisgroup.github.io/geneSetSimplifyR/reference/stopWords.md)
  : Stop Words
