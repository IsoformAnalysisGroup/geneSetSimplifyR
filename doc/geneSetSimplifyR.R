## ----setup, include=FALSE-----------------------------------------------------
rmdformats::downcute(default_style='dark')
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = TRUE
)

## ----installation, eval=FALSE-------------------------------------------------
#  if (!requireNamespace("devtools", quietly = TRUE)){
#    install.packages("devtools")
#  }
#  devtools::install_github("kvittingseerup/geneSetSimplifyR", build_vignettes = TRUE)

## ----load, message=FALSE------------------------------------------------------
library(geneSetSimplifyR)

data("exampleGeneSets")
data("exampleEnrichment")

## ----explore example data-----------------------------------------------------
### The gene sets are a list
class(exampleGeneSets)

### Where each entry is a gene set and its annoated genes
str(head(exampleGeneSets,2))

### The optional data have both enrichment_score and source
head( exampleEnrichment, 2 )

## ----geneSetSimplifyeR--------------------------------------------------------
exampleGsList <- geneSetSimplifyR(
  geneSetsList = exampleGeneSets,
  geneSetsDF = exampleEnrichment,
  resolution = c(0.1, 0.6, 1.1, 1.6), # for quick run - delete when running on real data
  verbose = FALSE
)

exampleGsList

## ----optimal cluster----------------------------------------------------------
extractResolution(exampleGsList)

## ----load example, message=FALSE----------------------------------------------
data("exampleGsList")
exampleGsList

## ----barchart, width="75%"----------------------------------------------------
plotBarchart(
  geneSetsList = exampleGsList,
  highligthClusterNo = c(0,1)
)

## ----umap highlight-----------------------------------------------------------
plotUMAP(
  geneSetsList = exampleGsList,
  highligthClusterNo = c(0,1),
  removeClusterId = TRUE,
  pointAlpha = 1 # because there are relatively few gene-sets here
)

## ----proportion, width="75%"--------------------------------------------------
plotProportions(
  geneSetsList = exampleGsList,
  highligthClusterNo = c(0,1),
  removeClusterId = TRUE
)

## ----enrichment, width="75%"--------------------------------------------------
plotEnrichment(
  geneSetsList = exampleGsList,
  highligthClusterNo = c(0,1),
  removeClusterId = TRUE
)

## ----violin, warning=FALSE, width="75%"---------------------------------------
plotViolin(
  geneSetsList = exampleGsList,
  highligthClusterNo = c(0,1),
  removeClusterId = FALSE
)

## ----dendogram, warning=FALSE, width="75%"------------------------------------
plotDendogram(
  geneSetsList = exampleGsList,
  removeClusterId = TRUE
)

## ----plotClustree-------------------------------------------------------------
plotClustree(exampleGsList)

## ----umap new res-------------------------------------------------------------
plotUMAP(
  geneSetsList = exampleGsList,
  resolution = 1.6,
  pointAlpha = 1 # because there are relatively few gene-sets here
)

## ----update res---------------------------------------------------------------
exampleGsList <- updateResolution(
  geneSetsList = exampleGsList,
  resolution = 1.6
)

exampleGsList

## ----plot labels, warning=FALSE-----------------------------------------------
plotClusterLabels(
  geneSetsList = exampleGsList,
  clusterNumber = 0,
  maxChar = 35 # shorten gene-set names to fit in vignette
) +
  theme_bw(base_size = 8) # make font smaller to fit in vignette

## ----update label, warning=FALSE----------------------------------------------
exampleGsList <- updateClusterLabels(
  geneSetsList = exampleGsList,
  clusterNumber = 0,
  newClusterName = 'cancer response'
)

plotUMAP(
  geneSetsList = exampleGsList,
  pointAlpha = 1 # because there are relatively few gene-sets here
)

## ----extract full table, warning=FALSE----------------------------------------
gsTable <- extractTable(
  geneSetsList = exampleGsList,
  summarize = FALSE
)

head(gsTable)

## ----extract summary table, warning=FALSE-------------------------------------
gsSummary <- extractTable(
  geneSetsList = exampleGsList,
  summarize = TRUE
)

head(gsSummary)

## ----ggplot2 example, width="75%"---------------------------------------------
data("exampleGsList")

p1 <- plotBarchart(
  geneSetsList = exampleGsList,
  highligthClusterNo = c(0,1)
)
class(p1)

p1

## ----ggplot2 mod, width="75%"-------------------------------------------------
p1 + 
  ggtitle('New plot tilte')

## ----ggplot2 mod2, width="75%"------------------------------------------------
p1 + 
  scale_fill_manual(values = c('black','green'))

## ----fgsea data---------------------------------------------------------------
library(fgsea)

data(examplePathways)
data(exampleRanks)

## ----example gene-sets--------------------------------------------------------
sigGenes <- tail(names(exampleRanks), 750)
head(sigGenes)

## ----fora---------------------------------------------------------------------
oraRes <- fora(
  pathways = examplePathways,
  genes    = sigGenes,
  universe = names(exampleRanks)
)
sum(oraRes$padj < 0.05)

## ----enrichment calculation---------------------------------------------------

n_genes <- length(sigGenes)
n_universe <- length(exampleRanks)


oraSig <-
  tibble::as_tibble(oraRes) %>% 
  dplyr::filter(
    padj < 0.05
  ) %>% 
  dplyr::mutate(
    enrichment_score = log2(
      (
        ( overlap / size       ) /
        ( n_genes / n_universe )
      ) + 0.06
    )
  )

nrow(oraSig)

head(oraSig)

## ----session info-------------------------------------------------------------
sessionInfo()

