#' Run LSI (truncated SVD) with a convergence fallback
#'
#' The workflow's LSI step normally runs through [Signac::RunSVD()]. On
#' difficult membership matrices the underlying solver can fail to converge and
#' throw an error ("did not converge -- try increasing work or maxit"). The
#' `irlba` backend Signac exposes (e.g. Signac 1.14.0) only lets us grow the
#' working subspace (`irlba.work`), and once that hits its ceiling (the smaller
#' matrix dimension) there is no lever left -- `maxit` is never plumbed through.
#'
#' This helper therefore first attempts the ordinary [Signac::RunSVD()] call and
#' returns its result unchanged when it succeeds. Only if that errors does it
#' fall back to computing the SVD directly with `RSpectra::svds()`, reproducing
#' Signac's LSI construction (feature loadings, singular-value derived `stdev`,
#' z-scored cell embeddings) so the resulting `lsi` reduction is a drop-in
#' replacement. `svds()` exposes both convergence levers through its `opts`
#' list: `ncv` (the Lanczos working subspace, analogous to `irlba.work`) and
#' `maxitr` (the iteration cap). The fallback retries with a progressively
#' larger `ncv` and `maxitr` until the decomposition converges.
#'
#' @param object Seurat object holding the (TF-IDF transformed) assay.
#' @param n Maximum number of singular vectors / LSI components to compute.
#'   Capped internally at one below the smaller matrix dimension.
#' @param reduction.name,reduction.key Name and key under which the reduction
#'   is stored (defaults match Signac's `"lsi"` / `"LSI_"`).
#' @param scale.max Optional clip value for the scaled embeddings.
#' @param verbose Logical; whether to message about progress and retries.
#'
#' @return The input Seurat object with the `lsi` reduction added.
#' @keywords internal
#' @noRd
runSVDsafe <- function(object,
                       n              = 50,
                       reduction.name = "lsi",
                       reduction.key  = "LSI_",
                       scale.max      = NULL,
                       verbose        = TRUE) {

  # Fast path: the ordinary Signac::RunSVD(). Return unchanged when it works.
  result <- tryCatch(
    suppressWarnings(
      Signac::RunSVD(object, approx = FALSE, verbose = FALSE)
    ),
    error = function(e) e
  )
  if (!inherits(result, "error")) {
    return(result)
  }

  if (verbose) {
    message("Signac::RunSVD() did not converge (", conditionMessage(result),
            "); falling back to RSpectra::svds() with escalating parameters...")
  }

  # Fallback path: compute the SVD directly via RSpectra::svds(), which exposes
  # both the working subspace ('ncv') and the iteration cap ('maxitr').
  runSVDviaSvds(
    object,
    n              = n,
    reduction.name = reduction.name,
    reduction.key  = reduction.key,
    scale.max      = scale.max,
    verbose        = verbose
  )
}

#' Compute LSI with RSpectra::svds(), retrying on non-convergence
#'
#' Reproduces [Signac::RunSVD()]'s LSI construction but drives the
#' decomposition with `RSpectra::svds()` so `ncv` and `maxitr` can be
#' escalated until it converges. See [runSVDsafe()] for context.
#'
#' @inheritParams runSVDsafe
#' @return The input Seurat object with the `lsi` reduction added.
#' @keywords internal
#' @noRd
runSVDviaSvds <- function(object,
                          n              = 50,
                          reduction.name = "lsi",
                          reduction.key  = "LSI_",
                          scale.max      = NULL,
                          verbose        = TRUE) {

  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("Package 'RSpectra' is required for the SVD convergence fallback. ",
         "Install it with install.packages('RSpectra').")
  }

  assay <- Seurat::DefaultAssay(object)

  # The matrix the SVD runs on: variable features x cells (gene sets), mirroring
  # Signac's RunSVD.Assay (subset to VariableFeatures, use the 'data' layer).
  features <- Seurat::VariableFeatures(object)
  if (length(features) == 0) {
    features <- rownames(object)
  }
  data.use <- Seurat::GetAssayData(object, assay = assay, layer = "data")[features, , drop = FALSE]

  # A = t(data.use): cells x features. svds needs k < ncv <= min(dim(A)).
  A      <- Matrix::t(data.use)
  maxNcv <- min(dim(A))                  # = min(nCells, nFeatures)
  n      <- min(n, ncol(data.use) - 1, maxNcv - 1)

  # Escalating (ncv, maxitr) attempts. ncv is svds' working subspace (like
  # irlba's 'work'); maxitr is the iteration cap the irlba path never exposed.
  # The first attempt stays close to svds' own defaults.
  ncvSteps    <- pmin(c(2 * n + 1, 3 * n, 5 * n, maxNcv), maxNcv)
  ncvSteps    <- pmax(ncvSteps, n + 1)   # ncv must exceed k = n
  maxitrSteps <- c(1000, 3000, 10000, 30000)

  components <- NULL
  lastError  <- NULL
  for (i in seq_along(ncvSteps)) {
    ncv    <- ncvSteps[i]
    maxitr <- maxitrSteps[i]

    if (verbose) {
      message(sprintf(
        "Running SVD via RSpectra::svds() with ncv = %d, maxitr = %d (attempt %d/%d)...",
        ncv, maxitr, i, length(ncvSteps)
      ))
    }

    components <- tryCatch(
      RSpectra::svds(
        A    = A,
        k    = n,
        opts = list(tol = 1e-5, ncv = ncv, maxitr = maxitr)
      ),
      error = function(e) {
        lastError <<- e
        NULL
      }
    )

    if (!is.null(components)) {
      break
    }
  }

  if (is.null(components)) {
    stop("SVD failed to converge after ", length(ncvSteps),
         " attempts with increasing 'ncv' and 'maxitr'. Last error: ",
         conditionMessage(lastError))
  }

  # Reproduce Signac's LSI construction (non-PCA path, scaled embeddings).
  feature.loadings <- components$v
  sdev             <- components$d / sqrt(max(1, nrow(data.use) - 1))
  cell.embeddings  <- components$u

  embed.mean      <- apply(cell.embeddings, MARGIN = 2, FUN = mean)
  embed.sd        <- apply(cell.embeddings, MARGIN = 2, FUN = stats::sd)
  norm.embeddings <- t((t(cell.embeddings) - embed.mean) / embed.sd)
  if (!is.null(scale.max)) {
    norm.embeddings[norm.embeddings >  scale.max] <-  scale.max
    norm.embeddings[norm.embeddings < -scale.max] <- -scale.max
  }

  rownames(feature.loadings) <- rownames(data.use)
  colnames(feature.loadings) <- paste0(reduction.key, seq_len(n))
  rownames(norm.embeddings)  <- colnames(data.use)
  colnames(norm.embeddings)  <- paste0(reduction.key, seq_len(n))

  reduction.data <- Seurat::CreateDimReducObject(
    embeddings = norm.embeddings,
    loadings   = feature.loadings,
    assay      = assay,
    stdev      = sdev,
    key        = reduction.key,
    misc       = components
  )

  object[[reduction.name]] <- reduction.data
  return(object)
}
