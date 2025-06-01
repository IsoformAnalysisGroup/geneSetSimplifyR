changeUsedResolution <- function(
    geneSetsList,
    resolution = NULL,
    removeClusterId = FALSE
) {
  if( ! is.null(resolution)) {
    resolutionColumn <-
      paste0("GSEA_snn_res.", resolution)
  } else {
    resolutionColumn <- geneSetsList@cluster@meta.data$best_cluster_res[1]
  }


  Seurat::Idents(geneSetsList@cluster) <-
    geneSetsList@cluster@meta.data[[resolutionColumn]]

  # get default labels
  defaultLabels <-
    geneSetsList@labels[[resolutionColumn]][["labelsDefault"]] %>%
    dplyr::mutate(cluster = as.numeric(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(weighted_score, n = 1, with_ties = F) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = paste0(cluster, ": ", label)) %>%
    dplyr::pull(label)

  if (removeClusterId) {
    defaultLabels <- stringr::str_remove(defaultLabels, "^[0-9]+: ")
  }

  currentClusterLabels <-
    levels(geneSetsList@cluster@active.ident)

  geneSetsList@cluster@active.ident <-
    plyr::mapvalues(x = geneSetsList@cluster@active.ident,
                    from = currentClusterLabels,
                    to = defaultLabels)
  return(geneSetsList)
}

#' Update default resolution of gsList
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#' @param resolution Float. The clustering resolution to use as the new default.
#' The new 'resolution' must have been used as a resolution in the 'clusterGeneSets()' function.
#'
#' @return 'gsList' object updated with new default resolution
#' @export
#'
#' @examples
#' gsListObject <- updateResolution(gsListObject)
updateResolution <- function(
    geneSetsList,
    resolution
) {
  geneSetsList <- changeUsedResolution(
    geneSetsList = geneSetsList,
    resolution = resolution,
    removeClusterId = FALSE
  )

  geneSetsList@cluster$best_cluster_res <-
    stringr::str_c('GSEA_snn_res.',resolution)

  return(geneSetsList)
}

#' Extract default resolution
#'
#' @param geneSetsList 'gsList' object that has been clustered with the 'clusterGeneSets()' function.
#'
#' @return default resolution
#' @export
#'
#' @examples
#' gsListObject <- updateResolution(gsListObject)
extractResolution <- function(
    geneSetsList
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

  geneSetsList@cluster$best_cluster_res[1] %>%
    stringr::str_remove('GSEA_snn_res.') %>%
    as.numeric() %>%
    return()
}


countConditions <- function(
    geneSetsList
) {
  nCond <- 'Not assigned'

  condCount <- dplyr::n_distinct(geneSetsList@enrichmentScore$source)

  if( condCount >= 3 ) {
    nCond <- '3+'
  }
  if( condCount == 1 ) {
    nCond <- '1'
  }
  if( condCount == 2 ) {
    nCond <- '2'
  }
  return(nCond)
}

cutGeneSetName <- function(aVec, maxChar = 60) {
  # aVec <- c(paste0(rep('a', 30), collapse = ''), paste0(rep('b', 50), collapse = ''), paste0(paste0(rep('c', 50), collapse = ''), '_up'))

  ### make vector with max lengths
  maxVecCar <- rep(maxChar, length(aVec)) - 3 # for the "..."

  ### Annotat up and down
  hasDirection <- which( grepl('_up$|_dn$|_down$', aVec ))
  if( length(hasDirection) ) {
    orgDirection <- sapply(
      strsplit(aVec[hasDirection], '_' ),
      function(x) tail(x,1)
    )

    # Update max lengths
    maxVecCar[hasDirection] <- maxVecCar[hasDirection] - (nchar(orgDirection) +1) # +1 for the "_"
  }


  ### Modify those which are to long
  toModify <- which(nchar(aVec) > maxVecCar)
  if(length(toModify)) {
    aVec[toModify] <- paste0(
      substr(aVec[toModify], 1, maxVecCar[toModify]),
      '...'
    )
  }

  ### Add up, down back
  addDirBack <- intersect(hasDirection, toModify)
  if( length(addDirBack) ) {
    aVec[addDirBack] <- paste0(
      aVec[addDirBack],
      '_',
      orgDirection
    )
  }

  return(aVec)
}

boldXaxisLabels <- function(src, levelsToBold) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- levelsToBold %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(levelsToBold, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'levelsToBold' must be in src")
  }
}
