## Shared test fixtures. Clustering is independent of the 'source' column, so
## objects with one or three sources are derived from the shipped example by
## editing its enrichment table.

smallInput <- function(n = 60) {
  sets <- exampleGeneSets[seq_len(n)]
  df <- as.data.frame(exampleEnrichment)
  df <- df[df$pathway %in% names(sets), , drop = FALSE]
  list(sets = sets, df = df)
}

oneSourceGsList <- function() {
  g <- exampleGsList
  es <- g@enrichmentScore
  es <- es[!duplicated(es$pathway), , drop = FALSE]
  es$source <- NA
  g@enrichmentScore <- es
  g
}

threeSourceGsList <- function() {
  g <- exampleGsList
  es <- g@enrichmentScore
  extra <- es[seq(1, nrow(es), by = 3), , drop = FALSE]
  extra$source <- "Third"
  extra$enrichment_score <- -extra$enrichment_score
  g@enrichmentScore <- rbind(es, extra)
  g
}

withoutPvalues <- function(g) {
  g@enrichmentScore$pval <- NA
  g@enrichmentScore$nlog10pval <- NA
  g
}

withPvalues <- function(g) {
  es <- g@enrichmentScore
  es$pval <- seq(1e-4, 0.05, length.out = nrow(es))
  es$nlog10pval <- -log10(es$pval)
  g@enrichmentScore <- es
  g
}

## Draw a plot to a null device so that errors in rendering surface
renderPlot <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  print(p)
  invisible(TRUE)
}
