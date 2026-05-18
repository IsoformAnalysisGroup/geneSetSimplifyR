#' @import ggplot2
#' @keywords internal
"_PACKAGE"

gsList <- setClass(
  "gsList",
  slots = c(
    genesetsDataframe = "data.frame",
    enrichmentScore = "data.frame",
    #source = "ANY",
    cluster = "ANY",
    labels = "list",
    metadata = "list"
  ),
  prototype = list(
    metadata = list()
  )
)
