### ---------------------------------------------------------------------------
### Pluggable clustering backends
###
### `clusterGeneSets()` always builds the same SNN graph (LSI -> cosine KNN ->
### shared-nearest-neighbor). Everything after that point -- turning the graph
### into cluster assignments -- is delegated to a *cluster method*, so the
### community-detection step can be swapped for another implementation (another
### package, a tuned variant, an in-house algorithm) without touching the rest
### of the workflow.
###
### A cluster method is a function with the contract
###
###     function(snn, resolution, verbose = FALSE, ...) -> cluster ids
###
###   snn        Shared-nearest-neighbour graph as a sparse symmetric matrix
###              (a Seurat 'Graph'); rows/columns are gene sets, values are SNN
###              similarities in [0, 1].
###   resolution A single numeric resolution.
###   verbose    Whether the backend may print progress.
###   ...        Extra arguments forwarded from `clusterGeneSets(clusterArgs=)`.
###
### The return value must be one cluster id per column of `snn`: a vector,
### factor, or single-column data.frame/matrix. Ids may be 0- or 1-based and in
### any order -- `clusterGeneSets()` renumbers them to start at 0. If the result
### is named, it is reordered to match `colnames(snn)`.
###
### Not every algorithm has a resolution knob. Methods that pick their own
### number of clusters (kmeans, apcluster, binary_cut) declare
### `usesResolution = FALSE`, and `clusterGeneSets()` then runs them once
### instead of sweeping a resolution vector that would produce identical
### results.
### ---------------------------------------------------------------------------

# Registry of available cluster methods, keyed by name. Populated in .onLoad().
.clusterMethodRegistry <- new.env(parent = emptyenv())


#' Define a clustering method for `clusterGeneSets()`
#'
#' @description Wraps a community-detection function into the object
#'   `clusterGeneSets()` expects. The function receives the shared-nearest-
#'   neighbour (SNN) graph built by `clusterGeneSets()` and a single resolution,
#'   and returns one cluster id per gene set.
#'
#' @param fn A function with signature `function(snn, resolution, verbose, ...)`
#'   returning one cluster id per column of `snn` (vector, factor, or
#'   single-column data.frame/matrix). See details.
#' @param name Short name for the method, used in messages and in the
#'   `algorithm` argument of `clusterGeneSets()`.
#' @param description One-line description of what the method does.
#' @param packages Character vector of packages the method needs. They are
#'   checked with `requireNamespace()` before the method runs, so a missing
#'   dependency produces a clear error instead of a cryptic one.
#' @param usesResolution Whether the method's output actually depends on
#'   `resolution`. Set to `FALSE` for algorithms that determine the number of
#'   clusters themselves (e.g. affinity propagation); `clusterGeneSets()` then
#'   runs the method once rather than repeating identical work for every
#'   requested resolution. Defaults to `TRUE`.
#'
#' @details `snn` is a sparse symmetric matrix (a Seurat `Graph`) whose
#'   rows/columns are gene sets and whose values are SNN similarities in
#'   `[0, 1]`. Returned ids may be 0- or 1-based; `clusterGeneSets()` renumbers
#'   them to start at 0 and orders the factor levels numerically. If the
#'   returned ids are named, they are reordered to match `colnames(snn)`.
#'
#' @return An object of class `gsClusterMethod`.
#' @export
#'
#' @seealso [registerClusterMethod()], [listClusterMethods()],
#'   [clusterGeneSets()]
#'
#' @examples
#' # A backend built on igraph's fast-greedy modularity optimization
#' fastGreedy <- clusterMethod(
#'   fn = function(snn, resolution, verbose = FALSE, ...) {
#'     g <- snnToIgraph(snn)
#'     igraph::membership(
#'       igraph::cluster_fast_greedy(g, weights = igraph::E(g)$weight)
#'     )
#'   },
#'   name = "fast_greedy",
#'   description = "igraph fast-greedy modularity (ignores resolution)",
#'   packages = "igraph",
#'   usesResolution = FALSE
#' )
#' fastGreedy
#'
#' # Try it on the SNN graph of the example data
#' ids <- fastGreedy$fn(exampleGsList@cluster@graphs$GSEA_snn, resolution = 1)
#' table(ids)
#'
#' # Use it in the workflow
#' geneSetsList <- initializeList(exampleGeneSets, exampleEnrichment,
#'                                verbose = FALSE)
#' gsListObject <- clusterGeneSets(geneSetsList, algorithm = fastGreedy,
#'                                 verbose = FALSE)
clusterMethod <- function(
    fn,
    name,
    description = NA_character_,
    packages = character(0),
    usesResolution = TRUE
) {
  if (!is.function(fn)) {
    stop("'fn' should be a function.")
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("'name' should be a single non-empty string.")
  }
  if (length(formals(fn)) < 2L && !"..." %in% names(formals(fn))) {
    stop(
      "'fn' should accept at least the arguments 'snn' and 'resolution' ",
      "(signature: function(snn, resolution, verbose = FALSE, ...))."
    )
  }

  structure(
    list(
      fn = fn,
      name = name,
      description = as.character(description)[1],
      packages = as.character(packages),
      usesResolution = isTRUE(usesResolution)
    ),
    class = "gsClusterMethod"
  )
}

#' @export
print.gsClusterMethod <- function(x, ...) {
  cat("<gsClusterMethod> ", x$name, "\n", sep = "")
  if (!is.na(x$description) && nzchar(x$description)) {
    cat("  ", x$description, "\n", sep = "")
  }
  if (length(x$packages)) {
    cat("  requires: ", paste(x$packages, collapse = ", "), "\n", sep = "")
  }
  if (!isTRUE(x$usesResolution)) {
    cat("  picks the number of clusters itself (resolution is ignored)\n")
  }
  invisible(x)
}


#' Register a clustering method for use by name
#'
#' @description Adds a clustering backend to the package registry so it can be
#'   selected with `clusterGeneSets(algorithm = "<name>")`. Registering is
#'   optional -- a `gsClusterMethod` (or a bare function) can also be passed
#'   straight to `clusterGeneSets()`.
#'
#' @param name Name to register the method under.
#' @param fn Either a `gsClusterMethod` (from [clusterMethod()]) or a function
#'   following the same contract; a bare function is wrapped automatically.
#' @param description One-line description; ignored when `fn` is already a
#'   `gsClusterMethod`.
#' @param packages Packages the method requires; ignored when `fn` is already a
#'   `gsClusterMethod`.
#' @param usesResolution Whether the method's output depends on `resolution`;
#'   ignored when `fn` is already a `gsClusterMethod`. See [clusterMethod()].
#' @param overwrite Whether to replace an already-registered method of the same
#'   name. Defaults to `FALSE`.
#'
#' @return The registered `gsClusterMethod`, invisibly.
#' @export
#'
#' @seealso [clusterMethod()], [listClusterMethods()]
#'
#' @examples
#' registerClusterMethod(
#'   "leiden_cpm",
#'   function(snn, resolution, verbose = FALSE, ...) {
#'     g <- snnToIgraph(snn)
#'     igraph::membership(igraph::cluster_leiden(
#'       g,
#'       objective_function = "CPM",
#'       resolution = resolution,
#'       weights = igraph::E(g)$weight
#'     ))
#'   },
#'   description = "igraph Leiden with the CPM objective",
#'   packages = "igraph",
#'   overwrite = TRUE
#' )
#' listClusterMethods()
#'
#' geneSetsList <- initializeList(exampleGeneSets, exampleEnrichment,
#'                                verbose = FALSE)
#' gsListObject <- clusterGeneSets(geneSetsList, algorithm = "leiden_cpm",
#'                                 resolution = 0.05, verbose = FALSE)
registerClusterMethod <- function(
    name,
    fn,
    description = NA_character_,
    packages = character(0),
    usesResolution = TRUE,
    overwrite = FALSE
) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("'name' should be a single non-empty string.")
  }
  if (!overwrite && exists(name, envir = .clusterMethodRegistry, inherits = FALSE)) {
    stop(
      "A clustering method named '", name, "' is already registered. ",
      "Use overwrite = TRUE to replace it."
    )
  }

  method <- if (inherits(fn, "gsClusterMethod")) {
    fn$name <- name
    fn
  } else {
    clusterMethod(
      fn = fn,
      name = name,
      description = description,
      packages = packages,
      usesResolution = usesResolution
    )
  }

  assign(name, method, envir = .clusterMethodRegistry)
  invisible(method)
}


#' List the available clustering methods
#'
#' @description Returns the clustering backends `clusterGeneSets()` can be
#'   pointed at by name, including any added with [registerClusterMethod()].
#'
#' @return A `data.frame` with columns `name`, `description`, `packages` and
#'   `uses_resolution` (`FALSE` for methods that choose the number of clusters
#'   themselves, which `clusterGeneSets()` therefore runs only once).
#' @export
#'
#' @seealso [clusterMethod()], [registerClusterMethod()]
#'
#' @examples
#' listClusterMethods()
listClusterMethods <- function() {
  names_ <- sort(ls(envir = .clusterMethodRegistry))
  methods_ <- lapply(names_, get, envir = .clusterMethodRegistry)

  data.frame(
    name = names_,
    description = vapply(methods_, function(m) m$description, character(1)),
    packages = vapply(
      methods_,
      function(m) paste(m$packages, collapse = ", "),
      character(1)
    ),
    uses_resolution = vapply(
      methods_,
      function(m) isTRUE(m$usesResolution),
      logical(1)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Resolve a clustering method specification
#'
#' @description Turns the `algorithm` argument of [clusterGeneSets()] -- a
#'   registered name, a `gsClusterMethod`, or a bare function -- into a
#'   `gsClusterMethod`.
#'
#' @param algorithm A registered method name, a `gsClusterMethod`, or a function
#'   following the cluster-method contract (see [clusterMethod()]).
#'
#' @return A `gsClusterMethod`.
#' @export
#'
#' @examples
#' getClusterMethod("louvain")
getClusterMethod <- function(algorithm) {
  if (inherits(algorithm, "gsClusterMethod")) {
    return(algorithm)
  }

  if (is.function(algorithm)) {
    return(clusterMethod(fn = algorithm, name = "custom"))
  }

  if (is.character(algorithm) && length(algorithm) == 1L) {
    if (exists(algorithm, envir = .clusterMethodRegistry, inherits = FALSE)) {
      return(get(algorithm, envir = .clusterMethodRegistry))
    }
    stop(
      "Unknown clustering method '", algorithm, "'. Valid options: ",
      paste(sort(ls(envir = .clusterMethodRegistry)), collapse = ", "),
      ". You can add your own with registerClusterMethod()."
    )
  }

  stop(
    "'algorithm' should be a method name, a gsClusterMethod object ",
    "(see clusterMethod()), or a function."
  )
}


#' Run a clustering method on an SNN graph
#'
#' Internal. Checks the method's dependencies, calls it for a single
#' resolution, and validates/orders the returned cluster ids.
#'
#' @param method A `gsClusterMethod`.
#' @param snn The SNN graph (sparse symmetric matrix).
#' @param resolution A single numeric resolution.
#' @param verbose Passed to the method.
#' @param clusterArgs Named list of extra arguments forwarded to the method.
#'
#' @return An unnamed vector of cluster ids, ordered as `colnames(snn)`.
#' @keywords internal
runClusterMethod <- function(
    method,
    snn,
    resolution,
    verbose = FALSE,
    clusterArgs = list()
) {
  missingPkgs <- method$packages[
    !vapply(method$packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missingPkgs)) {
    stop(
      "Clustering method '", method$name, "' requires the package(s): ",
      paste(missingPkgs, collapse = ", "), ". Please install them."
    )
  }

  # 'verbose' is only passed to methods that can take it, so a method with
  # the minimal signature function(snn, resolution) also works.
  fnArgs <- names(formals(method$fn))
  callArgs <- list(snn = snn, resolution = resolution)
  if (any(c("verbose", "...") %in% fnArgs)) {
    callArgs$verbose <- verbose
  }

  ids <- do.call(method$fn, c(callArgs, clusterArgs))

  # Accept the shapes community-detection functions commonly return.
  if (is.data.frame(ids) || is.matrix(ids)) {
    if (ncol(ids) != 1L) {
      stop(
        "Clustering method '", method$name, "' returned ", ncol(ids),
        " columns; expected one cluster id per gene set."
      )
    }
    rn <- rownames(ids)
    ids <- if (is.data.frame(ids)) ids[[1]] else ids[, 1]
    if (!is.null(rn)) names(ids) <- rn
  }
  if (inherits(ids, "communities")) {
    ids <- igraph::membership(ids)
  }

  if (is.null(ids) || length(ids) != ncol(snn)) {
    stop(
      "Clustering method '", method$name, "' returned ",
      if (is.null(ids)) "NULL" else paste(length(ids), "cluster ids"),
      "; expected ", ncol(snn), " (one per gene set)."
    )
  }

  # Reorder by name when the backend names its output after the gene sets.
  # Names that only partly match the gene sets point to a bug in the backend,
  # so they are rejected rather than silently assigned by position.
  idNames <- names(ids)
  if (!is.null(idNames) && any(idNames %in% colnames(snn))) {
    if (anyDuplicated(idNames) || !setequal(idNames, colnames(snn))) {
      stop(
        "Clustering method '", method$name, "' returned names that do not ",
        "match the gene sets one-to-one."
      )
    }
    ids <- ids[colnames(snn)]
  }

  ids <- as.vector(ids)
  if (is.list(ids)) {
    stop(
      "Clustering method '", method$name, "' returned a list; expected one ",
      "cluster id per gene set."
    )
  }
  if (anyNA(ids)) {
    stop("Clustering method '", method$name, "' returned NA cluster ids.")
  }

  return(unname(ids))
}


#' Convert an SNN graph to an igraph object
#'
#' @description Helper for writing igraph-based clustering methods. Turns the
#'   sparse symmetric SNN matrix handed to a cluster method into a weighted
#'   undirected `igraph` graph whose vertex names are the gene sets.
#'
#' @param snn The SNN graph passed to a cluster method (a sparse symmetric
#'   matrix / Seurat `Graph`).
#'
#' @return A weighted undirected `igraph` object.
#' @export
#'
#' @seealso [clusterMethod()]
#'
#' @examples
#' g <- snnToIgraph(exampleGsList@cluster@graphs$GSEA_snn)
#' g
snnToIgraph <- function(snn) {
  adjacency <- methods::as(snn, "dgCMatrix")

  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adjacency,
    mode = "max",        # symmetrise defensively; SNN is already symmetric
    weighted = TRUE,
    diag = FALSE
  )

  if (is.null(igraph::V(g)$name)) {
    igraph::V(g)$name <- colnames(snn)
  }

  return(g)
}


#' Convert an SNN graph to a dense similarity matrix
#'
#' @description Helper for writing clustering methods that operate on a full
#'   similarity matrix rather than a graph (affinity propagation, binary cut,
#'   k-means on similarity profiles). The SNN graph is already a similarity
#'   matrix with values in `[0, 1]`; this simply densifies it and keeps the gene
#'   set names as dimnames.
#'
#' @param snn The SNN graph passed to a cluster method.
#' @param warnAbove Number of gene sets above which a warning about the memory
#'   cost of densifying is emitted. Set to `Inf` to silence. The dense matrix
#'   needs roughly `8 * N^2` bytes, so 10,000 gene sets is about 800 MB.
#'
#' @return A dense symmetric numeric matrix with gene sets as rows/columns.
#' @export
#'
#' @seealso [clusterMethod()], [snnToIgraph()]
#'
#' @examples
#' simMat <- snnToMatrix(exampleGsList@cluster@graphs$GSEA_snn)
#' simMat[1:3, 1:3]
snnToMatrix <- function(snn, warnAbove = 5000) {
  n <- ncol(snn)
  if (n > warnAbove) {
    warning(
      "Densifying the SNN graph for ", n, " gene sets needs about ",
      round(8 * n^2 / 1024^3, 1),
      " GB. Consider a graph-based clustering method instead."
    )
  }

  simMatrix <- as.matrix(snn)
  dimnames(simMatrix) <- list(colnames(snn), colnames(snn))

  return(simMatrix)
}


### --- Built-in backends -----------------------------------------------------

# Seurat's FindClusters has a method for graphs, so the backend can hand it the
# SNN matrix directly rather than the whole Seurat object.
seuratBackend <- function(algo) {
  function(snn, resolution, verbose = FALSE, ...) {
    res <- tryCatch(
      Seurat::FindClusters(
        object = snn,
        resolution = resolution,
        algorithm = algo,
        verbose = verbose,
        ...
      ),
      error = function(e) {
        if (algo == 4) {
          stop(
            conditionMessage(e), " Seurat's Leiden algorithm may need extra ",
            "software (the Python package 'leidenalg' for older Seurat ",
            "versions). The 'igraph_leiden' method needs no extra software.",
            call. = FALSE
          )
        }
        stop(e)
      }
    )
    return(res)
  }
}

# igraph's community detection is randomized; 'seed' (default 42, NA for no
# seeding) makes it reproducible without changing the user's RNG stream.
igraphBackend <- function(fun, defaults = list()) {
  function(snn, resolution, verbose = FALSE, seed = 42, ...) {
    g <- snnToIgraph(snn)

    args <- utils::modifyList(
      defaults,
      list(weights = igraph::E(g)$weight, ...)
    )

    # igraph renamed the resolution argument across versions
    # (resolution_parameter -> resolution); support both.
    argNames <- names(formals(fun))
    if ("resolution" %in% argNames) {
      args$resolution <- resolution
    } else if ("resolution_parameter" %in% argNames) {
      args$resolution_parameter <- resolution
    } else {
      stop("The igraph clustering function has no resolution argument.")
    }

    withSeed(seed, igraph::membership(do.call(fun, c(list(g), args))))
  }
}

# Evaluate an expression under a fixed RNG seed without leaving the user's
# random stream perturbed. k-means initialization is random, so a seed is the
# difference between reproducible and non-reproducible cluster assignments.
withSeed <- function(seed, expr) {
  if (is.null(seed) || is.na(seed)) {
    return(force(expr))
  }
  withr::with_seed(seed, expr)
}

# Kneedle-style knee detection on the within-cluster-sum-of-squares curve,
# matching simplifyEnrichment::cluster_by_kmeans(). `elbowFinder` measures the
# perpendicular distance to the line joining the highest-WSS point and the
# largest-k point; `kneeFinder` measures the vertical gap below the chord
# spanning the whole curve. cluster_by_kmeans() takes the smaller of the two,
# which biases the choice toward fewer clusters.
elbowFinder <- function(x, y) {
  maxXx <- max(x)
  maxXy <- y[which.max(x)]
  maxYy <- max(y)
  maxYx <- x[which.max(y)]

  # A least-squares fit through two points is just the line through them.
  if (maxXx == maxYx) {
    return(x[1])
  }
  slope <- (maxXy - maxYy) / (maxXx - maxYx)
  intercept <- maxYy - slope * maxYx

  distances <- abs(slope * x - y + intercept) / sqrt(slope^2 + 1)
  return(x[which.max(distances)])
}

kneeFinder <- function(x, y) {
  n <- length(x)
  if (n < 2L || x[n] == x[1]) {
    return(x[1])
  }
  slope <- (y[n] - y[1]) / (x[n] - x[1])
  return(x[which.max(slope * x - y)])
}

# k-means over the rows of the SNN similarity matrix, i.e. gene sets are
# compared by their similarity profiles. k is chosen automatically the way
# simplifyEnrichment does it: sweep k, then take the more conservative of the
# elbow and knee estimates on the WSS curve.
kmeansBackend <- function(
    snn,
    resolution,
    verbose = FALSE,
    maxK = NULL,
    iterMax = 50,
    seed = 42,
    ...
) {
  simMatrix <- snnToMatrix(snn)
  n <- nrow(simMatrix)

  if (is.null(maxK)) {
    maxK <- max(2, min(round(n / 5), 100))
  }
  maxK <- min(maxK, n - 1L)
  if (maxK <= 2) {
    stop(
      "'maxK' should be larger than 2; got ", maxK,
      " (too few gene sets for automatic k selection)."
    )
  }

  ks <- 2:maxK
  clusterings <- vector("list", length(ks))
  wss <- numeric(length(ks))

  withSeed(seed, {
    for (i in seq_along(ks)) {
      km <- suppressWarnings(
        stats::kmeans(simMatrix, centers = ks[i], iter.max = iterMax, ...)
      )
      clusterings[[i]] <- km$cluster
      wss[i] <- sum(km$withinss)
    }
  })

  bestK <- min(elbowFinder(ks, wss), kneeFinder(ks, wss))

  if (verbose) {
    message("kmeans: selected k = ", bestK, " (searched k = 2 to ", maxK, ")")
  }

  ids <- clusterings[[which(ks == bestK)]]
  names(ids) <- rownames(simMatrix)
  return(ids)
}

# Affinity propagation. Unlike simplifyEnrichment, which turns the input into
# points and derives similarities with negDistMat(), the SNN graph already IS a
# similarity matrix, so it is handed to apcluster() directly -- keeping the
# "always cluster on SNN similarity" promise of the workflow.
#
# Granularity note: apcluster's usual `q` knob (the quantile of the input
# similarities used as the self-preference) is degenerate on an SNN graph.
# Roughly 80% of SNN entries are exactly 0, so every q below that quantile
# yields the same preference and the same clustering. Tune `p` instead --
# `clusterArgs = list(p = -1)` -- where lower (more negative) values relative to
# the [0, 1] similarity scale give fewer, larger clusters. Note that affinity
# propagation is invariant to adding a constant to all similarities and the
# preference together, so rescaling the matrix on its own changes nothing.
apclusterBackend <- function(
    snn,
    resolution,
    verbose = FALSE,
    ...
) {
  simMatrix <- snnToMatrix(snn)

  apRes <- apcluster::apcluster(
    s = simMatrix,
    details = FALSE,
    ...
  )

  if (!length(apRes@clusters)) {
    stop("apcluster() did not converge on any exemplars.")
  }

  ids <- integer(nrow(simMatrix))
  for (i in seq_along(apRes@clusters)) {
    ids[apRes@clusters[[i]]] <- i
  }
  names(ids) <- rownames(simMatrix)

  if (verbose) {
    message("apcluster: found ", length(apRes@clusters), " clusters")
  }

  return(ids)
}

# simplifyEnrichment's binary cut: recursively split the similarity matrix in
# two and stop when a node's similarity profile looks like a single coherent
# group. Operates on the similarity matrix, which the SNN graph already is.
binaryCutBackend <- function(
    snn,
    resolution,
    verbose = FALSE,
    cutoff = 0.85,
    ...
) {
  simMatrix <- snnToMatrix(snn)

  ids <- simplifyEnrichment::binary_cut(
    simMatrix,
    cutoff = cutoff,
    ...
  )
  names(ids) <- rownames(simMatrix)

  if (verbose) {
    message("binary_cut: found ", length(unique(ids)), " clusters")
  }

  return(ids)
}

.onLoad <- function(libname, pkgname) {
  register <- function(name, fn, description, packages, usesResolution = TRUE) {
    assign(
      name,
      clusterMethod(
        fn = fn,
        name = name,
        description = description,
        packages = packages,
        usesResolution = usesResolution
      ),
      envir = .clusterMethodRegistry
    )
  }

  register(
    "louvain", seuratBackend(1),
    "Seurat::FindClusters, original Louvain (algorithm 1)", "Seurat"
  )
  register(
    "louvain.multilevel", seuratBackend(2),
    "Seurat::FindClusters, Louvain with multilevel refinement (algorithm 2)",
    "Seurat"
  )
  register(
    "slm", seuratBackend(3),
    "Seurat::FindClusters, smart local moving (algorithm 3)", "Seurat"
  )
  register(
    "leiden", seuratBackend(4),
    "Seurat::FindClusters, Leiden via leidenalg (algorithm 4)",
    c("Seurat", "leiden")
  )
  register(
    "igraph_louvain",
    igraphBackend(igraph::cluster_louvain),
    "igraph::cluster_louvain on the SNN graph", "igraph"
  )
  register(
    "igraph_leiden",
    igraphBackend(
      igraph::cluster_leiden,
      defaults = list(objective_function = "modularity", n_iterations = 5)
    ),
    "igraph::cluster_leiden on the SNN graph (modularity objective)", "igraph"
  )

  ### Methods that determine the number of clusters themselves. They operate on
  ### the dense SNN similarity matrix rather than the graph, and ignore
  ### 'resolution' -- clusterGeneSets() runs them once.
  register(
    "kmeans", kmeansBackend,
    paste(
      "k-means on SNN similarity profiles; k chosen by elbow/knee on the WSS",
      "curve (as in simplifyEnrichment)"
    ),
    packages = character(0), usesResolution = FALSE
  )
  register(
    "apcluster", apclusterBackend,
    "Affinity propagation (apcluster::apcluster) on the SNN similarity matrix",
    packages = "apcluster", usesResolution = FALSE
  )
  register(
    "binary_cut", binaryCutBackend,
    "simplifyEnrichment::binary_cut on the SNN similarity matrix",
    packages = "simplifyEnrichment", usesResolution = FALSE
  )

  invisible()
}
