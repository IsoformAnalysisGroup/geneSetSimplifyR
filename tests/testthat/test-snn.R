snn <- exampleGsList@cluster@graphs$GSEA_snn

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

test_that("getSNN returns the SNN graph", {
    expect_identical(
        getSNN(exampleGsList), exampleGsList@cluster@graphs$GSEA_snn
    )
    unclustered <- initializeList(exampleGeneSets, verbose = FALSE)
    expect_error(getSNN(unclustered), "clusterGeneSets")
})
