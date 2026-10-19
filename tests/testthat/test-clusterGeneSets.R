test_that("the workflow reproduces the shipped example object", {
  g <- geneSetSimplifyR(
    exampleGeneSets, exampleEnrichment,
    removeFirstWord = TRUE, verbose = FALSE
  )
  resCols <- grep("^GSEA_snn_res", names(exampleGsList@cluster@meta.data),
                  value = TRUE)
  expect_equal(
    g@cluster@meta.data[, resCols],
    exampleGsList@cluster@meta.data[, resCols]
  )
  expect_equal(extractResolution(g), extractResolution(exampleGsList))
  expect_equal(
    levels(g@cluster@active.ident),
    levels(exampleGsList@cluster@active.ident)
  )
})

test_that("small inputs (50 or fewer gene sets) can be clustered", {
  input <- smallInput(30)
  g <- initializeList(input$sets, input$df, verbose = FALSE)
  g <- clusterGeneSets(g, resolution = c(0.6, 1.1), verbose = FALSE)
  expect_s4_class(g, "gsList")
  expect_lte(g@metadata$lsi_elbow_cutoff, 29)
  g <- labelClusters(g, verbose = FALSE)
  expect_s3_class(plotUMAP(g), "ggplot")
})

test_that("too few gene sets give a clear error", {
  input <- smallInput(8)
  g <- initializeList(input$sets, input$df, verbose = FALSE)
  expect_error(clusterGeneSets(g, verbose = FALSE), "At least 10 gene sets")
})

test_that("arguments are validated", {
  input <- smallInput(30)
  g <- initializeList(input$sets, input$df, verbose = FALSE)
  expect_error(clusterGeneSets(g, K = 1, verbose = FALSE), "'K'")
  expect_error(clusterGeneSets(g, K = 30, verbose = FALSE), "'K'")
  expect_error(clusterGeneSets(g, K = NULL, verbose = FALSE), "'K'")
  expect_error(clusterGeneSets(g, resolution = -1, verbose = FALSE), "positive")
  expect_error(
    clusterGeneSets(g, clusterArgs = list(resolution = 1), verbose = FALSE),
    "cannot contain"
  )
  expect_error(
    clusterGeneSets(g, clusterArgs = list(1), verbose = FALSE),
    "named"
  )
  expect_error(clusterGeneSets(g, elbow_method = "foo", verbose = FALSE))
})

test_that("clustering does not change the session's random numbers", {
  input <- smallInput(30)
  g <- initializeList(input$sets, input$df, verbose = FALSE)
  set.seed(1)
  expected <- stats::runif(3)
  set.seed(1)
  invisible(clusterGeneSets(g, resolution = 1, verbose = FALSE))
  expect_equal(stats::runif(3), expected)
})

test_that("an error lists why every resolution failed", {
  input <- smallInput(30)
  g <- initializeList(input$sets, input$df, verbose = FALSE)
  failing <- function(snn, resolution) stop("backend exploded")
  expect_error(
    clusterGeneSets(g, algorithm = failing, resolution = c(1, 2),
                    verbose = FALSE),
    "backend exploded"
  )
})

test_that("re-clustering removes stale labels", {
  input <- smallInput(30)
  g <- initializeList(input$sets, input$df, verbose = FALSE)
  g <- labelClusters(clusterGeneSets(g, resolution = 1, verbose = FALSE),
                     verbose = FALSE)
  expect_true(length(g@labels) > 0)
  g <- clusterGeneSets(g, resolution = 2, verbose = FALSE)
  expect_length(g@labels, 0)
})

test_that("the silhouette score handles degenerate clusterings", {
  d <- as.matrix(stats::dist(matrix(c(1, 2, 10, 11, 20), ncol = 1)))
  expect_equal(silhouetteScore(c(0, 0, 0, 0, 0), d), 0)
  expect_true(is.na(silhouetteScore(0:4, d)))
  expect_gt(silhouetteScore(c(0, 0, 1, 1, 2), d), 0)
})

test_that("the elbow works on short and flat scree curves", {
  expect_equal(kvsElbow(c(5, 4, 3))$elbow_cutoff, 3)
  expect_gte(kvsElbow(c(10, 5, 3, 2.5, 2.4, 2.3, 2.2))$elbow_cutoff, 2)
  expect_equal(kvsElbow(c(10, rep(1, 10)))$elbow_cutoff, 11)
  expect_gte(
    kvsElbow(c(10, rep(1, 20)), method = "zscore")$elbow_cutoff, 2
  )
  expect_error(kvsElbow(c(1, 2, 3)), "decreasing")
})
