#' Determine the elbow cutoff to be used to run UMAP
#'
#' @param variance Standard deviations of the LSI components, stored in a
#'   Seurat object as `seuratObj@reductions$lsi@stdev`. Must be sorted in
#'   decreasing order.
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
#' @details Short scree curves (fewer than 50 components, i.e. inputs with 50
#'   or fewer gene sets) are supported. With fewer than four components after
#'   dropping the first one, all components are used.
#'
#' @return A list with element `elbow_cutoff`: the index of the last LSI
#'   component to use (at least 2).
#' @keywords internal
kvsElbow <- function(
    variance,
    method = c("distance", "zscore"),
    zScoreCutoff = 2.5
) {
  method <- match.arg(method)

  ### Test input
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }
  if (length(variance) < 2) {
    stop("'variance' should hold at least 2 components.")
  }

  ### Too few components to locate an elbow: use them all
  if (length(variance) < 5) {
    return(list(elbow_cutoff = length(variance)))
  }

  if (method == "zscore") {
    ### Scale
    varianceScaled <- scale(variance, center = TRUE, scale = TRUE)

    ### Calculate stepwise change (1st order derivative)
    stepwiseDecrease <-
      varianceScaled[-length(varianceScaled)] -
      varianceScaled[-1]

    ### Calculate z-score from stepwise change
    stepwiseMad <- mad(stepwiseDecrease)
    if (!is.finite(stepwiseMad) || stepwiseMad == 0) {
      stepwiseZ <- rep(0, length(stepwiseDecrease))
    } else {
      stepwiseZ <-
        (stepwiseDecrease - median(stepwiseDecrease)) / stepwiseMad
    }

    if(!any(stepwiseZ > zScoreCutoff)) {
      cutoffIndex <- min(13L, length(variance))
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
    yRange <- max(y) - min(y)
    # A flat curve has no elbow (smoothing can leave numerical noise)
    if (yRange <= sqrt(.Machine$double.eps) * max(abs(y), 1)) {
      return(list(elbow_cutoff = length(variance)))
    }
    yNorm <- (y - min(y)) / yRange

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
