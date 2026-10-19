test_that("parseFgsea keeps significant pathways and p-values", {
  pathways <- names(exampleGeneSets)[1:30]
  res <- data.frame(
    pathway = c(pathways, "NOT_IN_SETS"),
    pval = c(rep(c(0.001, 0.5), 15), 0.001),
    padj = c(rep(c(0.01, 0.6), 15), 0.01),
    NES = seq(-2, 2, length.out = 31)
  )
  expect_warning(
    parsed <- parseFgsea(res, exampleGeneSets, cutoff = 0.05, source = "A"),
    "not in 'geneSets'"
  )
  expect_equal(nrow(parsed$geneSetsDF), 15)
  expect_true(all(c("pval", "padj") %in% names(parsed$geneSetsDF)))
  expect_named(parsed$geneSetsList, parsed$geneSetsDF$pathway)

  g <- initializeList(parsed$geneSetsList, parsed$geneSetsDF, verbose = FALSE)
  expect_false(all(is.na(g@enrichmentScore$nlog10pval)))
})

test_that("parseFgsea accepts ORA results and reports missing columns", {
  pathways <- names(exampleGeneSets)[1:10]
  ora <- data.frame(pathway = pathways, padj = 0.01, foldEnrichment = 2)
  parsed <- parseFgsea(ora, exampleGeneSets)
  expect_equal(parsed$geneSetsDF$enrichment_score[1], log2(2 + 1e-6))

  expect_error(
    parseFgsea(data.frame(pathway = pathways), exampleGeneSets),
    "missing required column"
  )
})

test_that("parsePairedGSEA makes one source per block", {
  pathways <- names(exampleGeneSets)[1:20]
  ora <- data.frame(
    pathway = pathways,
    padj_expression = rep(c(0.01, 0.5), 10),
    enrichment_score_expression = 1,
    padj_splicing = 0.01,
    enrichment_score_splicing = -1
  )
  parsed <- parsePairedGSEA(ora, exampleGeneSets, cutoff_paired = NULL)
  expect_equal(sort(unique(parsed$geneSetsDF$source)),
               c("expression", "splicing"))
  expect_equal(nrow(parsed$geneSetsDF), 30)
  expect_error(parsePairedGSEA(ora, exampleGeneSets), "padj_paired")
})
