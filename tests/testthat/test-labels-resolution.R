test_that("renaming twice keeps the original labels", {
  g <- updateClusterLabels(exampleGsList, clusterNumber = 0,
                           newClusterName = "First")
  g <- updateClusterLabels(g, clusterNumber = 1, newClusterName = "Second")

  original <- extractTable(exampleGsList, summarize = FALSE)
  renamed <- extractTable(g, summarize = FALSE)

  originalOf <- function(table, label) {
    unique(table$initial_cluster_label[table$cluster == label])
  }
  firstOriginal <- unique(original$cluster[renamed$cluster == "First"])
  secondOriginal <- unique(original$cluster[renamed$cluster == "Second"])
  expect_equal(originalOf(renamed, "First"), firstOriginal)
  expect_equal(originalOf(renamed, "Second"), secondOriginal)
})

test_that("renaming validates its input", {
  expect_error(
    updateClusterLabels(exampleGsList, clusterNumber = 99,
                        newClusterName = "x"),
    "does not exist"
  )
  expect_error(
    updateClusterLabels(exampleGsList, clusterNumber = 0, newClusterName = ""),
    "non-empty"
  )
  expect_error(
    updateClusterLabels(exampleGsList, resolution = 0.2, clusterNumber = 0,
                        newClusterName = "x"),
    "Available resolution"
  )
})

test_that("a rename with a resolution makes it the active resolution", {
  g <- updateClusterLabels(exampleGsList, resolution = 0.6, clusterNumber = 0,
                           newClusterName = "Renamed")
  expect_equal(extractResolution(g), 0.6)

  nClusters06 <- nlevels(exampleGsList@cluster@meta.data$GSEA_snn_res.0.6)
  expect_equal(nlevels(g@cluster@active.ident), nClusters06)

  # Plots and tables with resolution = NULL follow the active resolution
  tableActive <- extractTable(g)
  expect_true("Renamed" %in% tableActive$cluster)
  expect_equal(length(unique(tableActive$cluster)), nClusters06)
  expect_true(any(grepl("Renamed", levels(plotBarchart(g)$data$cluster))))

  # ... while an explicit resolution still wins
  nClusters11 <- nlevels(exampleGsList@cluster@meta.data$GSEA_snn_res.1.1)
  expect_equal(
    length(unique(extractTable(g, resolution = 1.1)$cluster)),
    nClusters11
  )

  # A later rename without resolution renames within the active resolution
  g <- updateClusterLabels(g, clusterNumber = 1, newClusterName = "Second")
  expect_true(all(c("Renamed", "Second") %in% extractTable(g)$cluster))
})

test_that("updateResolution changes the default", {
  g <- updateResolution(exampleGsList, 1.1)
  expect_equal(extractResolution(g), 1.1)
  expect_equal(g@cluster@meta.data$best_cluster_res[1], "GSEA_snn_res.1.1")
  expect_error(updateResolution(exampleGsList, 7), "not used")
})

test_that("clusters with the same label stay separate without ids", {
  labels <- levels(exampleGsList@cluster@active.ident)
  firstLabel <- sub("^0: ", "", labels[1])
  expect_warning(
    g <- updateClusterLabels(exampleGsList, clusterNumber = 1,
                             newClusterName = firstLabel),
    "already labeled"
  )
  table <- extractTable(g, removeClusterId = TRUE)
  expect_equal(
    length(unique(table$cluster)),
    nlevels(exampleGsList@cluster@active.ident)
  )
})

test_that("automatic labels are unique per resolution", {
  for (res in names(exampleGsList@labels)) {
    top <- exampleGsList@labels[[res]]$labelsDefault
    top <- top[order(top$cluster, -top$weighted_score), ]
    top <- top[!duplicated(top$cluster), ]
    expect_false(anyDuplicated(top$label) > 0)
  }
})

test_that("stop words are matched literally and may be empty", {
  g <- labelClusters(exampleGsList, stopWords = character(0), verbose = FALSE)
  expect_true(all(lengths(lapply(g@labels, `[[`, "labelsDefault")) > 0))
  g <- labelClusters(exampleGsList, stopWords = c("c++", "(", "cell"),
                     verbose = FALSE)
  topLabels <- levels(g@cluster@active.ident)
  expect_false(any(grepl("\\bCell\\b", topLabels)))
})

test_that("long gene-set names are shortened with their direction kept", {
  long <- paste0(strrep("A", 70), "_DN")
  short <- cutGeneSetName(c(long, "SHORT_UP", strrep("B", 50)), maxChar = 50)
  expect_equal(nchar(short[1]), 50)
  expect_match(short[1], "\\.\\.\\._DN$")
  expect_equal(short[2], "SHORT_UP")
  expect_equal(short[3], strrep("B", 50))
})

test_that("show prints a summary", {
  expect_output(show(exampleGsList), "gsList object with 197 gene sets")
})

test_that("objects without a metadata slot are upgraded", {
  g <- upgradeGsList(exampleGsList)
  expect_identical(g, exampleGsList)
})
