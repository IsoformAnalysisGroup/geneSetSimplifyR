#' Scatter plot of UMAP dimensionality reduction
#'
#' @description Graphs the output of the UMAP dimensionality reduction as a scatterplot
#' using the 'DimPlot()' function from the Seurat package.
#' Gene sets are colored and labelled by their cluster.
#' When there is only one source in the data (default),
#' each point will be shaped by the direction of the enrichment score of the gene set.
#' When there are two sources in the data,
#' each point will be shaped by the source indicating whether the gene set was found
#' significant in either or both sources.
#'
#' @param geneSetsList 'gsList' object with clusters and labels.
#' @param resolution The cluster resolution to plot. Default to NULL which means the optimal clustering found by clusterGeneSets().
#' @param removeClusterId Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.
#' @param highligthClusterNo A vector with the cluster numbers of the cluster (at the resolution) to highligth in the umap
#' @param highligthPointAlpha The transparency of the highligthed clusters. Set different from "alpha" to use transparency for highligting clusters. Default is 2/3
#' @param highligthPointSize The size of the highligthed clusters. Set different from "size" to use size for highligting clusters. Default is 2
#' @param highligthLabelAlpha The transparency of the highligthed clusters. Set different from "alpha" to use transparency for highligting clusters. Default is 2/3
#' @param highligthLabelSize The size of the highligthed clusters. Set different from "size" to use size for highligting clusters. Default is 2
#' @param pointSize size of points. Default is 2
#' @param pointAlpha transparency of points. Default is 1/3.
#' @param labelSize size of label Default is 2
#' @param labelAlpha transparency of label Default is 1.
#' @param interactive Logical, indicating whether to plot an interactive UMAP plot. Default: FALSE.
#' If TRUE, an interactive plot will be plotted using the 'HoverLocator()' from the Seurat package.
#'
#' @return A list of ggplot objects.
#' @export
#'
#' @examples plotUMAP(exampleGsList)
plotUMAP <- function(
    ### High level
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE,
    interactive = FALSE,
    highligthClusterNo = NULL,
    ### Detailed
    pointSize = 2,
    pointAlpha = 1/3,
    labelSize = 3,
    labelAlpha = 1,
    highligthPointAlpha = 2/3,
    highligthPointSize = 4,
    highligthLabelAlpha = 1,
    highligthLabelSize = 5
) {
  # Check input
  if (!class(geneSetsList) == "gsList") {
    stop("'geneSetsList' should be a gsList object")
  }

  if (length(geneSetsList@cluster) == 0) {
    stop(
      "No clustering performed for 'geneSetsList'.
         Please run 'clusterGeneSets()' and 'labelClusters()' before plotting."
    )
  }

  if (length(geneSetsList@labels) == 0) {
    stop(
      "No labels created for 'geneSetsList'.
         Please run 'labelClusters()' before plotting."
    )
  }

  hasEnrichmentScore <- ! all(is.na( geneSetsList@enrichmentScore$enrichment_score  ) )

  ### Extract clusters to highligth
  if(! is.null(highligthClusterNo)) {

    if( ! is.null(resolution)) {
      geneSetsList <- changeUsedResolution(
        geneSetsList = geneSetsList,
        resolution = resolution,
        removeClusterId = FALSE
      )
    }

    ### Extact cluster to highligth
    clusterPattern <- stringr::str_c('^',highligthClusterNo, ': ', collapse = '|')

    clusterNameToHighligth <-
      stringr::str_subset(
        string = levels(geneSetsList@cluster@active.ident),
        pattern = clusterPattern
      )

    if(removeClusterId) {
      clusterNameToHighligth <-
        stringr::str_remove(clusterNameToHighligth, "^[0-9]+: ")
    }
  }

  ### Update cluster resulution and remove ids
  if( ! is.null(resolution) | removeClusterId ) {
    geneSetsList <- changeUsedResolution(
      geneSetsList = geneSetsList,
      resolution = resolution,
      removeClusterId = removeClusterId
    )
  }

  ### Assign conditions
  nCond <- countConditions(geneSetsList)

  ### Many conditions
  if( nCond == '3+' ) {

    if (   interactive ) {
      pUMAP <- Seurat::HoverLocator(
        plot = pUMAP,
        information = Seurat::FetchData(
          geneSetsList@cluster,
          vars = c("ident",
                   "nCount_GSEA")
        )
      )
    }
    if ( ! interactive ) {
      pUMAP <- Seurat::DimPlot(
        geneSetsList@cluster,
        label = TRUE,
        label.box = TRUE,
        repel = TRUE,
        label.size = labelSize,
        pt.size = pointSize,
        group.by = 'ident'
      ) +
        theme_bw() +
        guides(colour = "none") +
        labs(title = NULL,x='UMAP1',y='UMAP2')
    }


  }

  ### Two conditions
  if( nCond == '2' ) {

    geneSetsDFFiltered <-
      geneSetsList@enrichmentScore %>%
      dplyr::count(pathway) %>%
      dplyr::right_join(
        x = .,
        y = geneSetsList@enrichmentScore[,c('source','pathway')],
        by = "pathway",
        multiple = "all"
      ) %>%
      dplyr::mutate(
        source = dplyr::case_when(
          n == 2 ~ "both",
          n == 1  ~ source
        )
      ) %>%
      dplyr::select(pathway, source) %>%
      dplyr::distinct()

    geneSetsList@cluster <-
      SeuratObject::AddMetaData(
        object = geneSetsList@cluster,
        metadata = geneSetsDFFiltered$source[match(
          geneSetsList@cluster@meta.data$gene_set, geneSetsDFFiltered$pathway
        )],
        col.name = 'source'
      )

    if (   interactive ) {
      pUMAP <- Seurat::HoverLocator(
        plot = pUMAP,
        information = Seurat::FetchData(
          geneSetsList@cluster,
          vars = c("ident",
                   "source",
                   "nCount_GSEA")
        )
      )
    }
    if ( ! interactive ) {
      pUMAP <- Seurat::DimPlot(
        geneSetsList@cluster,
        label = TRUE,
        label.box = TRUE,
        repel = TRUE,
        label.size = labelSize,
        pt.size = pointSize,
        shape.by = 'source',
        group.by = 'ident'
      ) +
        theme_bw() +
        theme(legend.position = "bottom") +
        guides(colour = "none") +
        labs(shape = "Gene set source", x='UMAP1',y='UMAP2') +
        labs(title = NULL)


    }



  }

  ### One condition
  if( nCond == '1' ) {

    # With enrichment scores
    if(   hasEnrichmentScore ) {
      enrichment_score_sign <-
        ifelse(geneSetsList@enrichmentScore$enrichment_score > 0, "Positive", "Negative")

      geneSetsList@cluster <-
        Seurat::AddMetaData(
          object = geneSetsList@cluster,
          metadata = enrichment_score_sign,
          col.name = 'enrichment_score_sign'
        )

      geneSetsList@cluster <-
        Seurat::AddMetaData(
          object = geneSetsList@cluster,
          metadata = geneSetsList@enrichmentScore$enrichment_score,
          col.name = 'enrichment_score'
        )

      if(   interactive ) {
        pUMAP <- Seurat::HoverLocator(
          plot = pUMAP,
          information = Seurat::FetchData(
            geneSetsList@cluster,
            vars = c("ident",
                     "enrichment_score",
                     "nCount_GSEA")
          )
        )
      }
      if( ! interactive ) {
        pUMAP <- Seurat::DimPlot(
          geneSetsList@cluster,
          label = TRUE,
          label.box = TRUE,
          repel = TRUE,
          label.size = labelSize,
          pt.size = pointSize,
          shape.by = 'enrichment_score_sign',
          group.by = 'ident'
        ) +
          theme_bw() +
          theme(legend.position = "bottom") +
          guides(colour = "none") +
          labs(shape = "Enrichment score") + #+ Seurat::NoLegend()
          labs(title = NULL, x='UMAP1',y='UMAP2')

      }
    }
    if( ! hasEnrichmentScore ) {
      if(   interactive ) {
        pUMAP <- Seurat::HoverLocator(
          plot = pUMAP,
          information = Seurat::FetchData(
            geneSetsList@cluster,
            vars = c("ident",
                     "nCount_GSEA")
          )
        )
      }
      if( ! interactive ) {
        pUMAP <- Seurat::DimPlot(
          geneSetsList@cluster,
          label = TRUE,
          label.box = TRUE,
          repel = TRUE,
          label.size = labelSize,
          pt.size = pointSize,
          group.by = 'ident'
        ) +
          theme_bw() +
          guides(colour = "none") +
          labs(title = NULL, x='UMAP1',y='UMAP2')

      }
    }

  }

  ### Set transparency and size
  if(! is.null(highligthClusterNo)) {

    ### Modify object
    # points
    pUMAP[[1]]$layers[[1]]$aes_params$alpha <-
      ifelse(
        test = pUMAP[[1]]$data$ident %in% clusterNameToHighligth,
        yes = highligthPointAlpha,
        no = pointAlpha
      )
    pUMAP[[1]]$layers[[1]]$aes_params$size <-
      ifelse(
        test = pUMAP[[1]]$data$ident %in% clusterNameToHighligth,
        yes = highligthPointSize,
        no = pointSize
      )

    # labels
    pUMAP[[1]]$layers[[2]]$aes_params$alpha <-
      ifelse(
        test = pUMAP[[1]]$layers[[2]]$data$ident %in% clusterNameToHighligth,
        yes = highligthLabelAlpha,
        no = labelAlpha
      )
    pUMAP[[1]]$layers[[2]]$aes_params$size <-
      ifelse(
        test = pUMAP[[1]]$layers[[2]]$data$ident %in% clusterNameToHighligth,
        yes =  highligthLabelSize,
        no =  labelSize
      )

  } else {
    pUMAP[[1]]$layers[[1]]$aes_params$alpha <- pointAlpha
    pUMAP[[1]]$layers[[1]]$aes_params$size <- pointSize

    # labels
    pUMAP[[1]]$layers[[2]]$aes_params$alpha <- labelAlpha
    pUMAP[[1]]$layers[[2]]$aes_params$size <- labelSize
  }

  ### Add title
  pUMAP <-
    pUMAP +
    labs(title = 'Gene Set Clustering')


  return(pUMAP)
}
