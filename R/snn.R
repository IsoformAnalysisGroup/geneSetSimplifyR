#' Get the SNN graph of a clustered gsList
#'
#' @description Returns the shared nearest neighbor (SNN) graph that
#'   `clusterGeneSets()` built and clustered. Convert it with [snnToIgraph()]
#'   or [snnToMatrix()] to analyze it with other tools.
#'
#' @param geneSetsList A 'gsList' object clustered with `clusterGeneSets()`.
#'
#' @return The SNN graph: a sparse symmetric matrix (a Seurat `Graph`) with the
#'   gene sets as rows and columns.
#' @export
#'
#' @seealso [snnToIgraph()], [snnToMatrix()]
#'
#' @examples
#' data(exampleGsList)
#' snn <- getSNN(exampleGsList)
#' dim(snn)
getSNN <- function(geneSetsList) {
    checkGsList(geneSetsList, needLabels = FALSE)
    geneSetsList@cluster@graphs$GSEA_snn
}


#' Convert an SNN graph to an igraph object
#'
#' @description Turns the SNN graph from [getSNN()] into a weighted
#'   undirected `igraph` graph whose vertex names are the gene sets.
#'
#' @param snn The SNN graph from [getSNN()] (a sparse symmetric matrix /
#'   Seurat `Graph`).
#'
#' @return A weighted undirected `igraph` object.
#' @export
#'
#' @seealso [getSNN()], [snnToMatrix()]
#'
#' @examples
#' data(exampleGsList)
#' g <- snnToIgraph(getSNN(exampleGsList))
#' g
snnToIgraph <- function(snn) {
    adjacency <- methods::as(snn, "dgCMatrix")

    g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = adjacency,
        mode = "max",        # symmetrize defensively; SNN is already symmetric
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
#' @description Turns the SNN graph from [getSNN()], whose values are
#'   similarities in `[0, 1]`, into a dense matrix with the gene-set names as
#'   dimnames, for methods that need a full similarity matrix.
#'
#' @param snn The SNN graph from [getSNN()].
#' @param warnAbove Number of gene sets above which a warning about the memory
#'   cost of densifying is emitted. Set to `Inf` to silence. The dense matrix
#'   needs roughly `8 * N^2` bytes, so 10,000 gene sets is about 800 MB.
#'
#' @return A dense symmetric numeric matrix with gene sets as rows/columns.
#' @export
#'
#' @seealso [getSNN()], [snnToIgraph()]
#'
#' @examples
#' data(exampleGsList)
#' simMat <- snnToMatrix(getSNN(exampleGsList))
#' simMat[1:3, 1:3]
snnToMatrix <- function(snn, warnAbove = 5000) {
    n <- ncol(snn)
    if (n > warnAbove) {
        warning(
            "Densifying the SNN graph for ", n, " gene sets needs about ",
            round(8 * n^2 / 1024^3, 1),
            " GB. Consider using the sparse graph instead."
        )
    }

    simMatrix <- as.matrix(snn)
    dimnames(simMatrix) <- list(colnames(snn), colnames(snn))

    return(simMatrix)
}
