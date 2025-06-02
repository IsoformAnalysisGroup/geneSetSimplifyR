#' Determine the elbow cutoff to be used to run UMAP
#'
#' @param variance Stored in a Seurat object as 'seuratObj@reductions$lsi@stde'
#' @param zScoreCutoff Integer indicating the cut off value for the Z-score
#'
#' @return Elbow cutoff
#'
#' @examples
#' # elbowRes <- kvsElbow(variance = seuratObj@reductions$lsi@stdev)
kvsElbow <- function(
    variance,
    zScoreCutoff=2.5
) {

  ### Test input
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }
  if(length(variance) < 50) {
    stop('variance should at least be 50 observations for this method to work')
  }

  ### Scale
  varianceScaled <- scale(variance, center = T, scale = TRUE)


  ### Calculate stepwise change (1st order derivative)
  stepwiseDecrease <-
    varianceScaled[-length(varianceScaled)] -
    varianceScaled[-1]


  ### Calculate z-score from stepwise change
  stepwiseZ <-
    (stepwiseDecrease - median(stepwiseDecrease)) /
    mad(stepwiseDecrease)


  if(!any(stepwiseZ > zScoreCutoff)) {
    cutoffIndex <- 13
  } else {
    ### Extract index of of cutoff
    cutoffIndex <- max(which(
      stepwiseZ > zScoreCutoff
    )) + 1
  }

  resList <- list(
    elbow_cutoff = cutoffIndex
  )

  return(resList)
}
