#' Plot dendrogram of cluster averages
#'
#' @description A plot showing the hierarchical clustering of the cluster
#'   averages (dendrogram). A tool made for showing how clusters (not
#'   individual gene-sets) relate to each other. The distance between two
#'   clusters is the average SNN distance (1 - SNN similarity) between their
#'   gene sets, and clusters are joined by complete linkage.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default `NULL` uses the
#'   active resolution: the optimal resolution found by `clusterGeneSets()`, or
#'   the one last set with `updateResolution()` or `updateClusterLabels()`.
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highlightClusterNo A vector with the cluster numbers of the clusters (at the resolution) to highlight (in bold) in the plot.
#'
#' @return A ggplot object. The underlying `hclust` object is stored in its
#'   `hclust` attribute.
#' @export
#'
#' @examples
#' plotDendrogram(exampleGsList)
#' hc <- attr(plotDendrogram(exampleGsList), "hclust")
plotDendrogram <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    highlightClusterNo = NULL
)  {
  # Check input
  checkGsList(geneSetsList)

  clusterNameToHighlight <- highlightedClusterNames(
    geneSetsList = geneSetsList,
    resolution = resolution,
    highlightClusterNo = highlightClusterNo,
    removeClusterId = removeClusterId
  )

  ### Set the cluster resolution and labels to plot
  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = resolution,
    removeClusterId = removeClusterId
  )

  ### Extract clusters
  clusterIdents <- geneSetsList@cluster@active.ident
  clusterOf <- as.character(clusterIdents)
  clusterVec <- levels(droplevels(clusterIdents))

  if (length(clusterVec) < 2) {
    stop("A dendrogram needs at least two clusters; this resolution has one.")
  }

  ### Average distance between the gene sets of each pair of clusters
  distMat <- as.matrix(1 - geneSetsList@cluster@graphs$GSEA_snn)
  distMat <- distMat[names(clusterIdents), names(clusterIdents)]

  avgDist <- matrix(
    0,
    nrow = length(clusterVec),
    ncol = length(clusterVec),
    dimnames = list(clusterVec, clusterVec)
  )
  for (i in seq_along(clusterVec)) {
    for (j in seq_along(clusterVec)) {
      if (i != j) {
        avgDist[i, j] <- mean(
          distMat[clusterOf == clusterVec[i], clusterOf == clusterVec[j]]
        )
      }
    }
  }

  # Perform hierarchical clustering using the complete linkage method
  hcAvg <- hclust(as.dist(avgDist), method = "complete")

  ### Draw the dendrogram as line segments, leaves on the y axis. The
  ### average distances between clusters are often all close to 1, so the
  ### axis starts just below the lowest merge instead of at 0.
  heightRange <- range(hcAvg$height)
  leafHeight <- heightRange[1] -
    max(diff(heightRange), 0.05 * heightRange[2], 1e-3)
  segments <- dendrogramSegments(hcAvg, leafHeight = leafHeight)
  leafLabels <- hcAvg$labels[hcAvg$order]

  p1 <-
    ggplot(segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(
      breaks = seq_along(leafLabels),
      labels = leafLabels,
      expand = expansion(add = 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(
      x = "Gene Set Cluster",
      y = "Average distance between clusters (1 - SNN)"
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  ### Make highlighted clusters bold
  if (length(clusterNameToHighlight)) {
    suppressWarnings(
      p1 <- p1 + theme(axis.text.y = element_text(
        face = ifelse(leafLabels %in% clusterNameToHighlight, "bold", "plain")
      ))
    )
  }

  attr(p1, "hclust") <- hcAvg
  return(p1)
}

#' Line segments of a dendrogram
#'
#' @param hc An `hclust` object.
#' @param leafHeight Height at which the leaves are drawn.
#' @return A data.frame with columns `x`, `y`, `xend` and `yend`; leaves are
#'   at `x = 1, 2, ...` in the order of `hc$order` and at `leafHeight`.
#' @keywords internal
dendrogramSegments <- function(hc, leafHeight = 0) {
  nLeaves <- length(hc$order)
  leafX <- integer(nLeaves)
  leafX[hc$order] <- seq_len(nLeaves)

  nodeX <- numeric(nLeaves - 1)
  nodeY <- hc$height
  position <- function(id) {
    if (id < 0) c(leafX[-id], leafHeight) else c(nodeX[id], nodeY[id])
  }

  segments <- vector("list", nLeaves - 1)
  for (i in seq_len(nLeaves - 1)) {
    a <- position(hc$merge[i, 1])
    b <- position(hc$merge[i, 2])
    h <- nodeY[i]
    nodeX[i] <- (a[1] + b[1]) / 2
    segments[[i]] <- data.frame(
      x = c(a[1], b[1], a[1]),
      y = c(a[2], b[2], h),
      xend = c(a[1], b[1], b[1]),
      yend = c(h, h, h)
    )
  }

  do.call(rbind, segments)
}
