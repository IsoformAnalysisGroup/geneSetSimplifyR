#' Cluster gene sets
#'
#' @param geneSetsList 'gsList' object created with the 'initializeList()'
#'   function.
#' @param resolution Vector with the desired clustering resolutions to use.
#'   Higher number usually results in more clusters.
#' @param K Integer or NA. Number of nearest neighbors used when building the
#'   KNN graph. If `NA` (default) it is set to `round(sqrt(N))`, where `N` is
#'   the number of gene sets. Must be smaller than `N`.
#' @param verbose a logic indicating whether to print progress statements.
#'   Defaults to TRUE.
#'
#' @details More than 50 gene sets are needed. The LSI uses up to 50
#'   dimensions (fewer when there are fewer genes); the first one,
#'   which mostly reflects gene-set size, is always left out. The gene sets
#'   are clustered with the Louvain algorithm (`Seurat::FindClusters()`) on
#'   the shared nearest neighbor (SNN) graph; see [getSNN()] to get the graph.
#'
#'   Clustering does not change the random number stream of the R session:
#'   the UMAP restores it afterwards.
#'   Re-clustering an object removes its existing labels, so run
#'   `labelClusters()` again afterwards.
#'
#' @return 'gsList' object with clustering stored. The best resolution (highest
#'  silhouette score) is annotated and used throughout the rest of the
#'  workflow (unless specified in subsequent functions).
#' @export
#'
#' @examples
#' data(exampleGeneSets, exampleEnrichment)
#' geneSetsList <- initializeList(
#'     geneSetsList = exampleGeneSets,
#'     geneSetsDF = exampleEnrichment,
#'     verbose = FALSE
#' )
#' gsListObject <- clusterGeneSets(
#'     geneSetsList = geneSetsList,
#'     resolution = c(0.6, 1.1),
#'     verbose = FALSE
#' )
clusterGeneSets <-
    function(
        geneSetsList,
        resolution = c(0.1, 0.3, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 4.1, 5.1),
        K = NA,
        verbose = TRUE
        ) {

        if (!inherits(geneSetsList, "gsList")) {
            stop("'geneSetsList' should be a gsList object.")
        }

        if (!is.numeric(resolution) || !length(resolution) ||
            anyNA(resolution) || any(!is.finite(resolution)) ||
            any(resolution <= 0)) {
            stop("'resolution' should be a vector of positive numbers.")
        }
        resolution <- unique(resolution)

        nGeneSets <- ncol(geneSetsList@genesetsDataframe)
        if (nGeneSets <= 50) {
            stop(
                "More than 50 gene sets are needed for clustering; got ",
                nGeneSets, "."
            )
        }

        if (length(K) != 1L) {
            stop("'K' should be NA or a single whole number.")
        }
        if (!is.na(K)) {
            if (!is.numeric(K) || K != round(K) || K < 2 || K >= nGeneSets) {
                stop(
                    "'K' should be NA or a whole number between 2 and the ",
                    "number of gene sets minus 1 (", nGeneSets - 1, ")."
                )
            }
        }

        if(verbose) {
            message("Clustering gene sets...")
        }

        suppressWarnings(
            seuratObj <-
                Seurat::CreateSeuratObject(
                    counts = geneSetsList@genesetsDataframe,
                    assay = "GSEA"
                )
        )

        # TF-IDF + truncated SVD; the TF-IDF method was chosen by grid search
        seuratObj <- runLSI(seuratObj, verbose = verbose)

        elbowRes <- kvsElbow(variance = seuratObj@reductions$lsi@stdev)

        # RunUMAP() reseeds the global RNG; keep the user's random numbers.
        # The option turns off the notice that older Seurat 5 releases give
        # about the default UMAP method; newer releases do not give it.
        seuratObj <- withr::with_preserve_seed(
            withr::with_options(
                list(Seurat.warn.umap.uwot = FALSE),
                Seurat::RunUMAP(
                    seuratObj,
                    reduction = "lsi",
                    dims = 2:elbowRes$elbow_cutoff,
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

        # dims, cosine metric and K were chosen by grid search
        seuratObj <- Seurat::FindNeighbors(
            object = seuratObj,
            reduction = 'lsi',
            dims = 2:elbowRes$elbow_cutoff,
            annoy.metric = "cosine",
            k.param = selectedK,
            compute.SNN = TRUE,
            verbose = FALSE
        )

        # Louvain (Seurat algorithm 1) on the SNN graph, one resolution at a
        # time, so a failing one (e.g. a very high one) is skipped
        snnGraph <- seuratObj@graphs$GSEA_snn

        resolutionsTested <- c()
        failures <- character(0)
        for (localRes in resolution) {
            clusters <- suppressWarnings(
                tryCatch(
                    Seurat::FindClusters(
                        object = snnGraph,
                        resolution = localRes,
                        algorithm = 1,
                        verbose = verbose
                    ),
                    error = function(e) e
                )
            )

            if (inherits(clusters, "error")) {
                failures[as.character(localRes)] <- conditionMessage(clusters)
                if (verbose) {
                    message(
                        "Resolution ", localRes, " failed (",
                        conditionMessage(clusters), ") - skipping."
                    )
                }
                next
            }

            # Ids start at 0 with numeric level order
            seuratObj@meta.data[[resolutionColumnName(localRes)]] <-
                numericClusterFactor(clusters[colnames(snnGraph), 1])

            resolutionsTested <- c(resolutionsTested, localRes)
        }

        if (!length(resolutionsTested)) {
            failureLines <- paste0(
                "  resolution ", names(failures), ": ", failures,
                collapse = "\n"
            )
            stop("Clustering failed for every resolution:\n", failureLines)
        }

        # Mean silhouette width on 1 - SNN; NA (undefined) is never chosen
        distMatrix <- as.matrix(1 - seuratObj@graphs$GSEA_snn)
        meanSil <- vapply(resolutionsTested, function(localRes) {
            localClust <- as.integer(as.character(
                seuratObj@meta.data[[resolutionColumnName(localRes)]]
            ))
            silhouetteScore(localClust, distMatrix)
        }, numeric(1))

        if (all(is.na(meanSil))) {
            bestRes <- resolutionsTested[1]
            warning(
                "The silhouette score could not be computed for any ",
                "resolution; using resolution ", bestRes, " as default.",
                call. = FALSE
            )
        } else {
            bestRes <- resolutionsTested[which.max(meanSil)]
        }

        Seurat::Idents(seuratObj) <-
            seuratObj@meta.data[[resolutionColumnName(bestRes)]]
        seuratObj@meta.data$best_cluster_res <- resolutionColumnName(bestRes)

        # Labels of an earlier clustering no longer match the new clusters
        geneSetsList@cluster <- seuratObj
        geneSetsList@labels <- list()

        silhouetteScores <-
            stats::setNames(meanSil, as.character(resolutionsTested))
        geneSetsList@metadata$k <- selectedK
        geneSetsList@metadata$lsi_elbow_cutoff <- elbowRes$elbow_cutoff
        geneSetsList@metadata$resolutions_tested <- resolutionsTested
        geneSetsList@metadata$selected_resolution <- bestRes
        geneSetsList@metadata$active_resolution <- resolutionColumnName(bestRes)
        geneSetsList@metadata$silhouette_scores <- silhouetteScores
        geneSetsList@metadata$clustered_at <- Sys.time()

        geneSetsList@cluster <-
            Seurat::AddMetaData(geneSetsList@cluster,
                geneSetsList@cluster@assays[["GSEA"]]$counts@Dimnames[[2]],
                col.name = "gene_set")

        if(verbose) {
            for (localRes in resolutionsTested) {
                resolutionColumn <- resolutionColumnName(localRes)
                clusterNumber <- length(levels(
                    geneSetsList@cluster@meta.data[[resolutionColumn]]
                ))

                message(
                    clusterNumber, " clusters found for resolution ", localRes
                )
            }
        }

        if(verbose) {
            message(
                "Determined ", bestRes,
                " as the optimal clustering resolution. This is set as default."
            )
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
