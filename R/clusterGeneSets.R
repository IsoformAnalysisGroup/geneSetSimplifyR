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
#' @param K Integer or NA. Number of nearest neighbours used when building the
#'   KNN graph. If `NA` (default) it is set to `min(round(sqrt(N)), 30)` where
#'   `N` is the number of gene sets.
#' @param elbow_method option: distance, zscore.
#' @param verbose a logic indicating whether to print progress statements. Defaults to TRUE.
#'
#' @return 'gsList' object with clustering stored. The best resolution (highest
#'  silhouette score) is annotated and used throughout the rest of the workflow (unless
#'  specified in subsequent functions).
#' @export
#'
#' @seealso [clusterMethod()], [registerClusterMethod()], [listClusterMethods()]
#'
#' @examples
#' # Create gsListObject
#' geneSetsList <- initializeList(geneSetsList = exampleGeneSets, geneSetsDF = exampleEnrichment)
#' gsListObject <- clusterGeneSets(geneSetsList = geneSetsList)
#'
#' # Swap in a different community-detection implementation
#' \dontrun{
#' gsListObject <- clusterGeneSets(geneSetsList, algorithm = "igraph_leiden")
#' }
clusterGeneSets <-
  function(
    geneSetsList,
    resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 3.1, 4.1, 5.1),
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

    if (!is.numeric(resolution)) {
      stop("'resolution' should be numeric.")
    }

    if (!is.list(clusterArgs)) {
      stop("'clusterArgs' should be a named list.")
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

    seuratObj <- Signac::RunTFIDF(seuratObj,
                                  verbose = FALSE,
                                  method = 1) # verified by new grid search

    seuratObj <- Signac::FindTopFeatures(seuratObj,
                                         min.cutoff = "q0")

    # LSI can fail to converge on difficult matrices; runSVDsafe() computes the
    # SVD via RSpectra::svds() and retries with larger 'ncv'/'maxitr' on failure.
    seuratObj <- runSVDsafe(seuratObj, verbose = verbose)

    # get max dim
    elbowRes <- kvsElbow(variance = seuratObj@reductions$lsi@stdev,
                          method = elbow_method)

    seuratObj <- suppressWarnings(
      Seurat::RunUMAP(
        seuratObj,
        reduction = "lsi",
        dims = 2:elbowRes$elbow_cutoff, # verified by new grid search
        verbose = FALSE
      )
    )
    
    assign_k <- function(N) {
      if(!is.na(K)){ return(K)}
      k <- round(sqrt(N))        # Take square root of N
      # if (k > 30) {       # Check if it exceeds 30
      #   k <- 30
      # }
      cat("Selected K:", k, "\n")
      return(k)
    }

    selectedK <- assign_k(nrow(seuratObj@meta.data))

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
          error = function(e) {
            if (verbose) {
              message(
                "Resolution ", localRes, " failed (", conditionMessage(e),
                ") - skipping."
              )
            }
            NULL
          }
        )
      )

      if (is.null(localIds)) {
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
      clusterResolution <- stringr::str_c('GSEA_snn_res.', localRes)
      seuratObj@meta.data[[clusterResolution]] <-
        numericClusterFactor(localIds)

      resolutionsTested <- c(resolutionsTested, localRes)
    }

    if (!length(resolutionsTested)) {
      message("No valid resolution. Please try with other resolutions.")
      return(NULL)
    }

    ### Set optimal clustering as default
    if(TRUE) {
      # Extract silluette scores
      meanSil <- sapply(resolutionsTested, function(localRes) {

        localClust <- as.integer(as.character(seuratObj@meta.data[,stringr::str_c('GSEA_snn_res.',localRes)]))

        if(length(unique(localClust)) > 1) {
          silInfo <- base::summary(
            cluster::silhouette(
              x = localClust,
              dmatrix = 1 - seuratObj@graphs$GSEA_snn
            )
          )
          mean_clust_sil <- mean(silInfo$clus.avg.widths)
          return(mean_clust_sil)
        } else {
          return(0)
        }

      })

      ### Extract best resolution
      bestRes <- resolutionsTested[which(meanSil == max(meanSil))[1] ]

      ### Set as default
      Seurat::Idents(seuratObj) <- seuratObj@meta.data[,stringr::str_c('GSEA_snn_res.',bestRes)]
      seuratObj@meta.data$best_cluster_res <- stringr::str_c('GSEA_snn_res.',bestRes)

    }

    # Save seurat object to list
    geneSetsList@cluster <- seuratObj

    ### Record clustering metadata
    silhouetteScores <- stats::setNames(meanSil, as.character(resolutionsTested))
    geneSetsList@metadata$k <- selectedK
    geneSetsList@metadata$lsi_elbow_cutoff <- elbowRes$elbow_cutoff
    geneSetsList@metadata$cluster_algorithm <- clusterMethodObj$name
    geneSetsList@metadata$resolutions_tested <- resolutionsTested
    geneSetsList@metadata$selected_resolution <- bestRes
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
