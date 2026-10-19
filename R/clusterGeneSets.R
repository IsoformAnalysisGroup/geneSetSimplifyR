#' Cluster gene sets
#'
#' @param geneSetsList 'gsList' object created with the 'initializeList()' function.
#' @param resolution Vector with the desired clustering resolutions to use. Higher number usually results in more clusters.
#' @param algorithm Clustering method applied to the SNN graph. Either the name
#'   of a registered method (see [listClusterMethods()]; e.g. 'louvain',
#'   'leiden', 'slm', 'igraph_louvain', 'igraph_leiden', 'kmeans', 'apcluster',
#'   'binary_cut'), a `gsClusterMethod` built with [clusterMethod()], or a
#'   function following the same contract. Only this step is swappable -- the
#'   SNN graph itself is always built the same way.
#'
#'   Methods that determine the number of clusters themselves ('kmeans',
#'   'apcluster', 'binary_cut') ignore `resolution` and are run once; the result
#'   is stored under the first requested resolution.
#' @param clusterArgs Named list of extra arguments passed on to the clustering
#'   method, for example `list(n.start = 20)` for the Seurat backends,
#'   `list(maxK = 50)` for 'kmeans', `list(p = -1)` for 'apcluster', or
#'   `list(cutoff = 0.9)` for 'binary_cut'.
#' @param K Integer or NA. Number of nearest neighbors used when building the
#'   KNN graph. If `NA` (default) it is set to `round(sqrt(N))`, where `N` is
#'   the number of gene sets. Must be smaller than `N`.
#' @param elbow_method Method for choosing the number of LSI dimensions from
#'   the scree curve: 'distance' (default) or 'zscore'. See the Details.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @details At least 10 gene sets are needed. The LSI uses up to 50
#'   dimensions (fewer when there are 50 gene sets or less); the first one,
#'   which mostly reflects gene-set size, is always left out.
#'
#'   Clustering does not change the random number stream of the R session:
#'   the UMAP and the seeded clustering backends restore it afterwards.
#'   Re-clustering an object removes its existing labels, so run
#'   `labelClusters()` again afterwards.
#'
#' @return 'gsList' object with clustering stored. The best resolution (highest
#'  silhouette score) is annotated and used throughout the rest of the workflow (unless
#'  specified in subsequent functions).
#' @export
#'
#' @seealso [clusterMethod()], [registerClusterMethod()], [listClusterMethods()]
#'
#' @examples
#' geneSetsList <- initializeList(
#'   geneSetsList = exampleGeneSets,
#'   geneSetsDF = exampleEnrichment,
#'   verbose = FALSE
#' )
#' gsListObject <- clusterGeneSets(
#'   geneSetsList = geneSetsList,
#'   resolution = c(0.6, 1.1),
#'   verbose = FALSE
#' )
#'
#' # Swap in a different community-detection implementation
#' gsListObject <- clusterGeneSets(
#'   geneSetsList,
#'   algorithm = "igraph_leiden",
#'   resolution = 1,
#'   verbose = FALSE
#' )
clusterGeneSets <-
  function(
    geneSetsList,
    resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 4.1, 5.1),
    algorithm = "louvain",
    clusterArgs = list(),
    K = NA,
    elbow_method = "distance",
    verbose = TRUE
) {

    # Check input
    if (!inherits(geneSetsList, "gsList")) {
      stop("'geneSetsList' should be a gsList object.")
    }
    geneSetsList <- upgradeGsList(geneSetsList)

    if (!is.numeric(resolution) || !length(resolution) ||
        anyNA(resolution) || any(!is.finite(resolution)) ||
        any(resolution <= 0)) {
      stop("'resolution' should be a vector of positive numbers.")
    }
    resolution <- unique(resolution)

    if (!is.list(clusterArgs)) {
      stop("'clusterArgs' should be a named list.")
    }
    if (length(clusterArgs) &&
        (is.null(names(clusterArgs)) || any(!nzchar(names(clusterArgs))))) {
      stop("Every element of 'clusterArgs' should be named.")
    }
    reservedArgs <- intersect(
      names(clusterArgs), c("snn", "resolution", "verbose")
    )
    if (length(reservedArgs)) {
      stop(
        "'clusterArgs' cannot contain ", quoteAndCollapse(reservedArgs),
        "; these are set by clusterGeneSets()."
      )
    }

    elbow_method <- match.arg(elbow_method, c("distance", "zscore"))

    nGeneSets <- ncol(geneSetsList@genesetsDataframe)
    if (nGeneSets < 10) {
      stop(
        "At least 10 gene sets are needed for clustering; got ", nGeneSets, "."
      )
    }

    if (length(K) != 1L) {
      stop("'K' should be NA or a single whole number.")
    }
    if (!is.na(K)) {
      if (!is.numeric(K) || K != round(K) || K < 2 || K >= nGeneSets) {
        stop(
          "'K' should be NA or a whole number between 2 and the number of ",
          "gene sets minus 1 (", nGeneSets - 1, ")."
        )
      }
    }

    # Resolve the clustering backend up front so a bad spec fails before the
    # expensive LSI/UMAP steps run.
    clusterMethodObj <- getClusterMethod(algorithm)

    # Some backends (kmeans, apcluster, binary_cut) pick the number of clusters
    # themselves. Sweeping resolutions there would just repeat identical work,
    # so collapse to a single run.
    if (!isTRUE(clusterMethodObj$usesResolution) && length(resolution) > 1) {
      if (verbose) {
        message(
          "'", clusterMethodObj$name, "' determines the number of clusters ",
          "itself; ignoring all but the first resolution (", resolution[1], ")."
        )
      }
      resolution <- resolution[1]
    }

    if(verbose) {
      message("Clustering gene sets...")
    }

    # Create seurat object
    suppressWarnings(
      seuratObj <-
        Seurat::CreateSeuratObject(counts = geneSetsList@genesetsDataframe,
                                   assay = "GSEA")
    )

    # LSI: TF-IDF (method 1, verified by grid search) followed by a truncated
    # SVD. runLSI() falls back to RSpectra::svds() if the SVD does not converge.
    seuratObj <- runLSI(seuratObj, verbose = verbose)

    # get max dim
    elbowRes <- kvsElbow(variance = seuratObj@reductions$lsi@stdev,
                          method = elbow_method)

    # Seurat::RunUMAP() seeds the global RNG itself (seed.use = 42); keep the
    # user's random number stream as it was.
    seuratObj <- withr::with_preserve_seed(
      suppressWarnings(
        Seurat::RunUMAP(
          seuratObj,
          reduction = "lsi",
          dims = 2:elbowRes$elbow_cutoff, # verified by new grid search
          verbose = FALSE
        )
      )
    )

    if (is.na(K)) {
      selectedK <- round(sqrt(nGeneSets))
    } else {
      selectedK <- as.integer(K)
    }
    if (verbose) {
      message("Using K = ", selectedK, " nearest neighbors.")
    }

    seuratObj <- Seurat::FindNeighbors(
      object = seuratObj,
      reduction = 'lsi',
      dims = 2:elbowRes$elbow_cutoff, # verified by new grid search
      annoy.metric = "cosine",        # verified by new grid search
      k.param = selectedK,            # updated by new grid search
      compute.SNN = TRUE,
      verbose = FALSE
    )

    ### Find clusters
    ### The SNN graph is the hand-off point to the clustering backend: every
    ### method sees the same graph and returns one cluster id per gene set.
    ### Resolutions are run one at a time so a single failing resolution is
    ### dropped instead of taking the whole run with it.
    snnGraph <- seuratObj@graphs$GSEA_snn

    resolutionsTested <- c()
    failures <- character(0)
    for (localRes in resolution) {
      localIds <- suppressWarnings(
        tryCatch(
          runClusterMethod(
            method = clusterMethodObj,
            snn = snnGraph,
            resolution = localRes,
            verbose = verbose,
            clusterArgs = clusterArgs
          ),
          error = function(e) e
        )
      )

      if (inherits(localIds, "error")) {
        failures[as.character(localRes)] <- conditionMessage(localIds)
        if (verbose) {
          message(
            "Resolution ", localRes, " failed (", conditionMessage(localIds),
            ") - skipping."
          )
        }
        next
      }

      ### Make sure the cluster number starts from 0

      ### Normalize every resolution's cluster factor once, at the source:
      ### renumber so ids start at 0 AND force numeric level order (0,1,2,...).
      ### Backends differ here -- Seurat is 0-based, igraph is 1-based, and
      ### Seurat sometimes stores the levels lexicographically ("10" sorts
      ### before "2"), which silently breaks any downstream code that pairs
      ### labels to factor levels by position. Fixing it here keeps every
      ### consumer correct without per-function guards.
      ### `normalizeClusterLevels()` is the same implementation, exported for
      ### migrating objects clustered earlier.
      clusterResolution <- resolutionColumnName(localRes)
      seuratObj@meta.data[[clusterResolution]] <-
        numericClusterFactor(localIds)

      resolutionsTested <- c(resolutionsTested, localRes)
    }

    if (!length(resolutionsTested)) {
      stop(
        "Clustering failed for every resolution with method '",
        clusterMethodObj$name, "':\n",
        paste0("  resolution ", names(failures), ": ", failures,
               collapse = "\n")
      )
    }

    ### Set optimal clustering as default
    # Mean silhouette width per resolution on the SNN distance (1 - SNN).
    # A single cluster scores 0. A resolution where the silhouette is not
    # defined (e.g. every gene set in its own cluster) scores NA and is never
    # chosen.
    distMatrix <- as.matrix(1 - seuratObj@graphs$GSEA_snn)
    meanSil <- vapply(resolutionsTested, function(localRes) {
      localClust <- as.integer(as.character(
        seuratObj@meta.data[[resolutionColumnName(localRes)]]
      ))
      silhouetteScore(localClust, distMatrix)
    }, numeric(1))

    ### Extract best resolution
    if (all(is.na(meanSil))) {
      bestRes <- resolutionsTested[1]
      warning(
        "The silhouette score could not be computed for any resolution; ",
        "using resolution ", bestRes, " as default.",
        call. = FALSE
      )
    } else {
      bestRes <- resolutionsTested[which.max(meanSil)]
    }

    ### Set as default
    Seurat::Idents(seuratObj) <-
      seuratObj@meta.data[[resolutionColumnName(bestRes)]]
    seuratObj@meta.data$best_cluster_res <- resolutionColumnName(bestRes)

    # Save seurat object to list. Labels of an earlier clustering no longer
    # match the new clusters, so they are removed.
    geneSetsList@cluster <- seuratObj
    geneSetsList@labels <- list()

    ### Record clustering metadata
    silhouetteScores <- stats::setNames(meanSil, as.character(resolutionsTested))
    geneSetsList@metadata$k <- selectedK
    geneSetsList@metadata$lsi_elbow_cutoff <- elbowRes$elbow_cutoff
    geneSetsList@metadata$cluster_algorithm <- clusterMethodObj$name
    geneSetsList@metadata$resolutions_tested <- resolutionsTested
    geneSetsList@metadata$selected_resolution <- bestRes
    geneSetsList@metadata$active_resolution <- resolutionColumnName(bestRes)
    geneSetsList@metadata$silhouette_scores <- silhouetteScores
    geneSetsList@metadata$clustered_at <- Sys.time()

    # Add meta data column with gene set names
    geneSetsList@cluster <-
      Seurat::AddMetaData(geneSetsList@cluster,
                          geneSetsList@cluster@assays[["GSEA"]]$counts@Dimnames[[2]],
                          col.name = "gene_set")

    if(verbose) {
      for (localRes in resolutionsTested) {
        resolutionColumn <- paste0("GSEA_snn_res.", localRes)
        clusterNumber <-
          length(levels(geneSetsList@cluster@meta.data[[resolutionColumn]]))

        message(paste0(clusterNumber, " clusters found for resolution ", localRes))
      }
    }

    if(verbose) {
      message(paste(
        'Determined',
        bestRes,
        "as the optimal clustering resolution. This is set as default.",
        sep = ' '
      ))
    }

    if(verbose) {
      message("Clustering DONE")
    }

    return(geneSetsList)
  }

#' Mean silhouette width of a clustering
#'
#' Internal helper for `clusterGeneSets()`.
#'
#' @param clusters Integer cluster ids, one per gene set.
#' @param distMatrix Dense distance matrix between the gene sets.
#' @return The mean of the per-cluster average silhouette widths; 0 for a
#'   single cluster and `NA` when the silhouette is not defined.
#' @keywords internal
silhouetteScore <- function(clusters, distMatrix) {
  nClusters <- length(unique(clusters))
  if (nClusters == 1) {
    return(0)
  }
  if (nClusters >= length(clusters)) {
    return(NA_real_)
  }

  sil <- tryCatch(
    cluster::silhouette(x = clusters, dmatrix = distMatrix),
    error = function(e) NA
  )
  if (!inherits(sil, "silhouette")) {
    return(NA_real_)
  }

  mean(summary(sil)$clus.avg.widths)
}
