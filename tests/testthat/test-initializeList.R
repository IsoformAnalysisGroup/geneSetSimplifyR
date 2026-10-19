test_that("gene sets without a table are initialized", {
  g <- initializeList(exampleGeneSets[1:20], verbose = FALSE)
  expect_s4_class(g, "gsList")
  expect_equal(ncol(g@genesetsDataframe), 20)
  expect_true(all(is.na(g@enrichmentScore$enrichment_score)))
  expect_true(all(is.na(g@enrichmentScore$nlog10pval)))
  expect_true(all(is.na(g@enrichmentScore$source)))
})

test_that("NES is used as enrichment score", {
  input <- smallInput(20)
  df <- input$df
  names(df)[names(df) == "enrichment_score"] <- "NES"
  g <- initializeList(input$sets, df, verbose = FALSE)
  expect_equal(g@enrichmentScore$enrichment_score, df$NES)
  expect_false("NES" %in% names(g@enrichmentScore))
})

test_that("pval and pvalue give the same -log10 p-values", {
  input <- smallInput(20)
  df <- input$df
  df$pval <- seq(0, 0.5, length.out = nrow(df))
  g1 <- initializeList(input$sets, df, verbose = FALSE)
  names(df)[names(df) == "pval"] <- "pvalue"
  g2 <- initializeList(input$sets, df, verbose = FALSE)
  expect_equal(g1@enrichmentScore$nlog10pval, g2@enrichmentScore$nlog10pval)
  expect_true(all(is.finite(g1@enrichmentScore$nlog10pval)))
  expect_equal(g1@enrichmentScore$nlog10pval[2], -log10(df$pvalue[2]))
})

test_that("invalid p-values are rejected", {
  input <- smallInput(20)
  df <- input$df
  df$pval <- 2
  expect_error(initializeList(input$sets, df, verbose = FALSE), "between 0 and 1")
  df$pval <- "small"
  expect_error(initializeList(input$sets, df, verbose = FALSE), "numeric")
})

test_that("duplicated pathway/source rows are rejected", {
  input <- smallInput(20)
  df <- rbind(input$df, input$df[1, ])
  df$enrichment_score[nrow(df)] <- 99
  expect_error(initializeList(input$sets, df, verbose = FALSE), "only once")

  # Exact duplicates are simply removed
  df <- rbind(input$df, input$df[1, ])
  g <- initializeList(input$sets, df, verbose = FALSE)
  expect_equal(nrow(g@enrichmentScore), nrow(input$df))
})

test_that("gene-set names are validated", {
  sets <- exampleGeneSets[1:5]
  expect_error(initializeList(unname(sets), verbose = FALSE), "named")
  dup <- c(sets, sets[1])
  expect_error(initializeList(dup, verbose = FALSE), "unique")
})

test_that("missing pathways are named in the error", {
  input <- smallInput(20)
  df <- input$df
  df$pathway[1] <- "NOT_A_GENE_SET"
  expect_error(initializeList(input$sets, df, verbose = FALSE), "NOT_A_GENE_SET")
})

test_that("empty gene sets are removed with a warning", {
  sets <- exampleGeneSets[1:12]
  sets$EMPTY_SET <- character(0)
  sets$ONLY_NA <- NA_character_
  expect_warning(
    g <- initializeList(sets, verbose = FALSE),
    "contain no genes"
  )
  expect_false(any(c("EMPTY_SET", "ONLY_NA") %in% colnames(g@genesetsDataframe)))
  expect_false(any(c("EMPTY_SET", "ONLY_NA") %in% g@enrichmentScore$pathway))
})

test_that("list columns such as leadingEdge are dropped", {
  input <- smallInput(20)
  df <- input$df
  df$leadingEdge <- I(lapply(seq_len(nrow(df)), function(i) c("A", "B")))
  g <- initializeList(input$sets, df, verbose = FALSE)
  expect_false("leadingEdge" %in% names(g@enrichmentScore))
})

test_that("reserved columns are removed with a warning", {
  input <- smallInput(20)
  df <- input$df
  df$cluster <- 1
  expect_warning(
    g <- initializeList(input$sets, df, verbose = FALSE),
    "reserved"
  )
  expect_false("cluster" %in% names(g@enrichmentScore))
})

test_that("gene sets can be read from a GMT file", {
  sets <- exampleGeneSets[1:15]
  gmtFile <- tempfile(fileext = ".gmt")
  writeLines(
    vapply(names(sets), function(name) {
      paste(c(name, "", sets[[name]]), collapse = "\t")
    }, character(1)),
    gmtFile
  )
  fromFile <- initializeList(gmtFile, verbose = FALSE)
  fromList <- initializeList(sets, verbose = FALSE)
  expect_equal(fromFile@genesetsDataframe, fromList@genesetsDataframe)
  expect_error(initializeList("no_such_file.gmt", verbose = FALSE),
               "File not found")
})
