snn <- exampleGsList@cluster@graphs$GSEA_snn

test_that("built-in methods are registered", {
  methods <- listClusterMethods()
  expect_true(all(c("louvain", "leiden", "igraph_louvain", "igraph_leiden",
                    "kmeans", "apcluster", "binary_cut") %in% methods$name))
  expect_false(methods$uses_resolution[methods$name == "kmeans"])
  expect_error(getClusterMethod("nope"), "Unknown clustering method")
})

test_that("igraph and kmeans backends return one id per gene set", {
  for (name in c("louvain", "igraph_louvain", "igraph_leiden", "kmeans")) {
    ids <- runClusterMethod(getClusterMethod(name), snn, resolution = 1)
    expect_length(ids, ncol(snn))
  }
})

test_that("seeded backends are reproducible", {
  m <- getClusterMethod("igraph_louvain")
  expect_equal(
    runClusterMethod(m, snn, resolution = 1),
    runClusterMethod(m, snn, resolution = 1)
  )
})

test_that("optional backends run when installed", {
  skip_if_not_installed("apcluster")
  ids <- runClusterMethod(getClusterMethod("apcluster"), snn, resolution = 1)
  expect_length(ids, ncol(snn))
})

test_that("a method with the minimal signature works", {
  minimal <- function(snn, resolution) rep(c(1, 2), length.out = ncol(snn))
  ids <- runClusterMethod(getClusterMethod(minimal), snn, resolution = 1)
  expect_length(ids, ncol(snn))
})

test_that("method output shapes are normalized", {
  oneColumn <- function(snn, resolution, ...) {
    matrix(rep(1:3, length.out = ncol(snn)), ncol = 1,
           dimnames = list(colnames(snn), NULL))
  }
  ids <- runClusterMethod(getClusterMethod(oneColumn), snn, resolution = 1)
  expect_equal(unique(ids), 1:3)

  # Named output is reordered to the gene sets
  reversed <- function(snn, resolution, ...) {
    ids <- seq_len(ncol(snn))
    names(ids) <- colnames(snn)
    rev(ids)
  }
  ids <- runClusterMethod(getClusterMethod(reversed), snn, resolution = 1)
  expect_equal(ids, seq_len(ncol(snn)))
})

test_that("malformed method output is rejected", {
  partialNames <- function(snn, resolution, ...) {
    ids <- rep(1, ncol(snn))
    names(ids) <- c(colnames(snn)[-1], "unknown")
    ids
  }
  expect_error(
    runClusterMethod(getClusterMethod(partialNames), snn, resolution = 1),
    "one-to-one"
  )
  tooShort <- function(snn, resolution, ...) 1:3
  expect_error(
    runClusterMethod(getClusterMethod(tooShort), snn, resolution = 1),
    "expected"
  )
  withNA <- function(snn, resolution, ...) c(NA, rep(1, ncol(snn) - 1))
  expect_error(
    runClusterMethod(getClusterMethod(withNA), snn, resolution = 1),
    "NA"
  )
})

test_that("cluster ids become consecutive and start at 0", {
  expect_equal(levels(numericClusterFactor(c(5, 1, 3, 3))), c("0", "1", "2"))
  expect_equal(
    as.character(numericClusterFactor(c(5, 1, 3, 3))),
    c("2", "0", "1", "1")
  )
  expect_equal(
    as.character(numericClusterFactor(c("b", "a", "b"))),
    c("1", "0", "1")
  )
  expect_equal(levels(numericClusterFactor(c(10, 2, 0))), c("0", "1", "2"))
})

test_that("snn helpers convert the graph", {
  expect_s3_class(snnToIgraph(snn), "igraph")
  m <- snnToMatrix(snn)
  expect_equal(dim(m), dim(snn))
  expect_equal(rownames(m), colnames(snn))
})
