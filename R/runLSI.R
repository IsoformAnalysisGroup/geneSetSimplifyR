#' Latent semantic indexing of the gene-set membership matrix
#'
#' Internal. Reproduces the LSI of `Signac` (`RunTFIDF(method = 1)`,
#' `FindTopFeatures(min.cutoff = "q0")` and `RunSVD()`) without depending on
#' Signac and its genomics dependencies:
#' \enumerate{
#'   \item TF-IDF: term frequency per gene set times inverse document frequency
#'     per gene, then `log1p(x * 1e4)`. Stored in the `data` layer.
#'   \item All genes are used, ordered by how many gene sets contain them.
#'   \item Truncated SVD with `irlba`, up to `n` components; embeddings are
#'     scaled per component and `stdev` is the singular value divided by
#'     `sqrt(number of genes - 1)`.
#' }
#' If `irlba` does not converge (it reports this as a warning), the SVD is
#' recomputed with `RSpectra::svds()` when RSpectra is installed, with
#' escalating `ncv`/`maxitr`, or else with base `svd()`.
#'
#' @param object Seurat object with the gene x gene-set counts.
#' @param n Maximum number of LSI components. Capped at the number of gene
#'   sets minus 1.
#' @param verbose Logical; whether to message about fallbacks.
#' @return The Seurat object with the TF-IDF `data` layer and an `lsi`
#'   reduction.
#' @keywords internal
#' @noRd
runLSI <- function(object, n = 50, verbose = TRUE) {
  assay <- Seurat::DefaultAssay(object)
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  counts <- methods::as(counts, "CsparseMatrix")

  ### TF-IDF (Signac method 1)
  normData <- tfidfTransform(counts)
  object <- SeuratObject::SetAssayData(
    object, assay = assay, layer = "data", new.data = normData
  )

  ### All features, ordered by total count (Signac min.cutoff = "q0")
  featureCounts <- Matrix::rowSums(counts)
  features <- names(featureCounts)[order(featureCounts, decreasing = TRUE)]
  SeuratObject::VariableFeatures(object, assay = assay) <- features

  dataUse <- normData[features, , drop = FALSE]
  n <- min(n, ncol(dataUse) - 1)

  components <- truncatedSVD(Matrix::t(dataUse), n = n, verbose = verbose)

  ### Embeddings and loadings as in Signac::RunSVD()
  reductionKey <- "LSI_"
  featureLoadings <- components$v
  sdev <- components$d / sqrt(max(1, nrow(dataUse) - 1))
  cellEmbeddings <- components$u

  embedMean <- apply(cellEmbeddings, MARGIN = 2, FUN = mean)
  embedSd <- apply(cellEmbeddings, MARGIN = 2, FUN = stats::sd)
  normEmbeddings <- t((t(cellEmbeddings) - embedMean) / embedSd)

  rownames(featureLoadings) <- rownames(dataUse)
  colnames(featureLoadings) <- paste0(reductionKey, seq_len(n))
  rownames(normEmbeddings) <- colnames(dataUse)
  colnames(normEmbeddings) <- paste0(reductionKey, seq_len(n))

  object[["lsi"]] <- Seurat::CreateDimReducObject(
    embeddings = normEmbeddings,
    loadings = featureLoadings,
    assay = assay,
    stdev = sdev,
    key = reductionKey,
    misc = components
  )

  return(object)
}

#' TF-IDF transform of a sparse count matrix (Signac method 1)
#'
#' @param counts Sparse genes x gene-sets matrix of 0/1 memberships.
#' @param scaleFactor Scale factor applied before `log1p()`.
#' @return Sparse matrix of the same shape.
#' @keywords internal
#' @noRd
tfidfTransform <- function(counts, scaleFactor = 1e4) {
  npeaks <- Matrix::colSums(counts)
  if (any(npeaks == 0)) {
    stop("Some gene sets contain no genes.")
  }
  tf <- Matrix::tcrossprod(x = counts, y = Matrix::Diagonal(x = 1 / npeaks))

  rsums <- Matrix::rowSums(counts)
  idf <- ncol(counts) / rsums

  normData <- Matrix::Diagonal(n = length(idf), x = idf) %*% tf
  normData@x <- log1p(normData@x * scaleFactor)
  normData@x[is.na(normData@x)] <- 0
  dimnames(normData) <- dimnames(counts)

  methods::as(normData, "CsparseMatrix")
}

#' Truncated SVD with a convergence fallback
#'
#' @param A Gene sets x genes matrix.
#' @param n Number of singular vectors.
#' @param verbose Logical; whether to message about fallbacks.
#' @return A list with `d`, `u` and `v`, as returned by `irlba::irlba()`.
#' @keywords internal
#' @noRd
truncatedSVD <- function(A, n, verbose = TRUE) {
  ### Small problems: a full SVD is exact and cheap. irlba is not designed
  ### for computing most of the singular values.
  if (2 * n >= min(dim(A))) {
    return(denseSVD(A, n))
  }

  ### irlba reports non-convergence as a warning ("did not converge"); treat
  ### that as a failure so the fallback runs, and drop other warnings.
  result <- tryCatch(
    withCallingHandlers(
      irlba::irlba(A = A, nv = n, work = n * 3, tol = 1e-05),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) {
          stop(conditionMessage(w), call. = FALSE)
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (!inherits(result, "error")) {
    return(result)
  }

  if (requireNamespace("RSpectra", quietly = TRUE)) {
    if (verbose) {
      message(
        "The LSI did not converge with irlba (", conditionMessage(result),
        "); retrying with RSpectra::svds()..."
      )
    }
    return(svdsWithRetries(A, n = n, verbose = verbose))
  }

  if (verbose) {
    message(
      "The LSI did not converge with irlba (", conditionMessage(result),
      "); computing a full SVD instead. Installing the 'RSpectra' package ",
      "makes this fallback faster for large inputs."
    )
  }
  denseSVD(A, n)
}

#' Exact SVD via base::svd(), truncated to n components
#' @keywords internal
#' @noRd
denseSVD <- function(A, n) {
  full <- svd(as.matrix(A), nu = n, nv = n)
  list(d = full$d[seq_len(n)], u = full$u, v = full$v)
}

#' Truncated SVD with RSpectra::svds(), retrying on non-convergence
#'
#' `svds()` exposes both convergence levers: `ncv` (the Lanczos working
#' subspace) and `maxitr` (the iteration cap). They are increased step by step
#' until the decomposition converges.
#' @keywords internal
#' @noRd
svdsWithRetries <- function(A, n, verbose = TRUE) {
  maxNcv <- min(dim(A))
  ncvSteps <- pmin(c(2 * n + 1, 3 * n, 5 * n, maxNcv), maxNcv)
  ncvSteps <- pmax(ncvSteps, n + 1) # ncv must exceed k = n
  maxitrSteps <- c(1000, 3000, 10000, 30000)

  lastError <- NULL
  for (i in seq_along(ncvSteps)) {
    if (verbose) {
      message(sprintf(
        "Running RSpectra::svds() with ncv = %d, maxitr = %d (attempt %d/%d)...",
        ncvSteps[i], maxitrSteps[i], i, length(ncvSteps)
      ))
    }

    components <- tryCatch(
      RSpectra::svds(
        A = A,
        k = n,
        opts = list(tol = 1e-5, ncv = ncvSteps[i], maxitr = maxitrSteps[i])
      ),
      error = function(e) e
    )

    if (!inherits(components, "error")) {
      return(components)
    }
    lastError <- components
  }

  stop(
    "The SVD failed to converge after ", length(ncvSteps), " attempts with ",
    "increasing 'ncv' and 'maxitr': ", conditionMessage(lastError)
  )
}
