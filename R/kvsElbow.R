#' Determine the elbow cutoff to be used to run UMAP
#'
#' @param variance Stored in a Seurat object as 'seuratObj@reductions$lsi@stde'
#' @param method Which elbow-detection method to use. One of:
#'   \describe{
#'     \item{`"distance"` (default)}{Returns the point of the (smoothed) scree
#'       curve that lies farthest from the chord connecting its endpoints. Light
#'       spline smoothing keeps it from being thrown off by jitter in the noisy
#'       tail, and it needs no threshold.}
#'     \item{`"zscore"`}{Flags the stepwise drops (1st derivative) that are
#'       statistical outliers via a robust z-score and returns the last such
#'       drop. Sensitive to the `zScoreCutoff` threshold.}
#'   }
#' @param zScoreCutoff Integer indicating the cut off value for the Z-score.
#'   Only used when `method = "zscore"`.
#'
#' @return Elbow cutoff
#'
#' @examples
#' # elbowRes <- kvsElbow(variance = seuratObj@reductions$lsi@stdev)
kvsElbow <- function(
    variance,
    method = c("distance", "zscore"),
    zScoreCutoff=2.5
) {
  method <- match.arg(method)

  ### Test input
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }
  if(length(variance) < 50) {
    stop('variance should at least be 50 observations for this method to work')
  }

  if (method == "zscore") {
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
      cutoffIndex <- max(max(which(stepwiseZ > zScoreCutoff)) + 1L, 2L)
    }
  } else {
    ### Distance-to-chord (Kneedle-style).
    ### LSI_1 typically captures technical variation (depth) and dominates the
    ### scree curve, so drop it before locating the knee; its inflated stdev
    ### would otherwise distort the [0, 1] normalization and the chord.
    varianceElbow <- variance[-1]

    n <- length(varianceElbow)
    x <- seq_len(n)

    ### Light smoothing to suppress point-wise jitter in the tail
    y <- stats::smooth.spline(x, varianceElbow)$y

    ### Normalize both axes to [0, 1] so neither dominates the distance
    xNorm <- (x - min(x)) / (max(x) - min(x))
    yNorm <- (y - min(y)) / (max(y) - min(y))

    ### Endpoints of the scree curve
    p1 <- c(xNorm[1], yNorm[1])
    pN <- c(xNorm[n], yNorm[n])

    ### Unit direction of the chord from first to last endpoint
    chord <- pN - p1
    chord <- chord / sqrt(sum(chord^2))

    ### Perpendicular distance of every point to that chord
    vecFromP1 <- cbind(xNorm - p1[1], yNorm - p1[2])
    proj <- as.vector(vecFromP1 %*% chord)
    perp <- vecFromP1 - outer(proj, chord)
    dist <- sqrt(rowSums(perp^2))

    ### +1 maps the index back to the original numbering (LSI_1 was dropped),
    ### so the downstream `dims = 2:elbow_cutoff` still lines up.
    cutoffIndex <- max(which.max(dist) + 1L, 2L)
  }

  resList <- list(
    elbow_cutoff = cutoffIndex
  )

  return(resList)
}
