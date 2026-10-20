twoSources <- exampleGsList
oneSource <- oneSourceGsList()
threeSources <- threeSourceGsList()
objects <- list(one = oneSource, two = twoSources, three = threeSources)

test_that("ggplot-based plots work for 1, 2 and 3+ sources", {
  for (name in names(objects)) {
    g <- objects[[name]]
    plots <- list(
      plotUMAP(g, highlightClusterNo = 0),
      plotUMAP(g, removeClusterId = TRUE, resolution = 0.6),
      plotBarchart(g, highlightClusterNo = 1),
      plotEnrichments(g, highlightClusterNo = 0, removeClusterId = TRUE),
      plotViolin(g, highlightClusterNo = 0),
      plotViolin(g, removeClusterId = TRUE, resolution = 1.1),
      plotDendrogram(g, highlightClusterNo = 0),
      plotClusterLabels(g, clusterNumber = 0),
      plotClustree(g, resolutionsToPlot = c(0.1, 0.6, 1.1))
    )
    for (p in plots) {
      expect_s3_class(p, "ggplot")
      expect_true(renderPlot(p))
    }
  }
})

test_that("plotProportions needs several sources", {
  expect_error(plotProportions(oneSource), "at least two sources")
  for (g in list(twoSources, threeSources)) {
    p <- plotProportions(g, highlightClusterNo = 0)
    expect_true(renderPlot(p))
    # Each gene set is counted once
    expect_equal(nrow(p$data), length(unique(g@enrichmentScore$pathway)))
  }
  expect_true("Multiple" %in% plotProportions(threeSources)$data$source)
})

test_that("plotHeatmap works for 1, 2 and 3+ sources", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (g in objects) {
    expect_s3_class(plotHeatmap(g, drawPlot = FALSE), "pheatmap")
    expect_s3_class(
      plotHeatmap(g, transpose = TRUE, highlightClusterNo = 0,
                  drawPlot = FALSE),
      "pheatmap"
    )
  }
  # Forced scaling of a single source does not fail
  expect_s3_class(plotHeatmap(oneSource, scaleCols = TRUE, drawPlot = FALSE),
                  "pheatmap")
  expect_error(plotHeatmap(twoSources, scaleCols = "yes"), "scaleCols")
})

test_that("plotHeatmap copes with constant values", {
  g <- twoSources
  g@enrichmentScore$enrichment_score <- 1
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(plotHeatmap(g, scaleCols = FALSE, drawPlot = FALSE),
                  "pheatmap")
})

test_that("plotFrom prefers -log10 p-values when present", {
  noPvalues <- withoutPvalues(twoSources)
  expect_equal(plotEnrichments(twoSources)$labels$y, "Mean -log10(p-value)")
  expect_equal(plotEnrichments(noPvalues)$labels$y, "Mean Enrichment")
  expect_equal(plotViolin(twoSources)$labels$y, "-log10(p-value)")
  expect_equal(plotViolin(noPvalues)$labels$y, "Enrichment")
  expect_equal(
    plotEnrichments(twoSources, plotFrom = "enrichment_score",
                    summaryStat = "median")$labels$y,
    "Median Enrichment"
  )
  expect_error(plotEnrichments(noPvalues, plotFrom = "nlog10pval"),
               "all values are NA")
  # Other numeric columns supplied by the user can be plotted too
  expect_equal(plotEnrichments(twoSources, plotFrom = "padj")$labels$y,
               "Mean padj")
})

test_that("the example data carries fgsea p-values", {
  expect_true(all(c("pval", "padj") %in% names(exampleEnrichment)))
  expect_true(all(exampleEnrichment$padj < 0.05))
  expect_equal(
    exampleGsList@enrichmentScore$nlog10pval,
    -log10(exampleGsList@enrichmentScore$pval)
  )
})

test_that("plotViolin keeps the median order without cluster ids", {
  withIds <- plotViolin(twoSources)
  withoutIds <- plotViolin(twoSources, removeClusterId = TRUE)
  expect_equal(
    sub("^[0-9]+: ", "", levels(withIds$data$clusterLabel)),
    levels(withoutIds$data$clusterLabel)
  )
})

test_that("plotDendrogram returns its clustering", {
  p <- plotDendrogram(twoSources)
  expect_s3_class(attr(p, "hclust"), "hclust")
})

test_that("interactive UMAP works", {
  skip_if_not_installed("plotly")
  p <- plotUMAP(twoSources, interactive = TRUE)
  expect_s3_class(p, "plotly")
})

test_that("plots give clear errors for unknown clusters and resolutions", {
  expect_error(plotUMAP(twoSources, resolution = 0.25), "Available resolution")
  expect_error(plotClusterLabels(twoSources, clusterNumber = 99),
               "clusterNumber")
  expect_warning(plotBarchart(twoSources, highlightClusterNo = 99),
                 "cannot be highlighted")
})

test_that("plotting a resolution does not change the object", {
  before <- twoSources@cluster@active.ident
  invisible(plotUMAP(twoSources, resolution = 1.1))
  expect_identical(twoSources@cluster@active.ident, before)
})
