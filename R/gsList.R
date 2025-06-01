gsList <- setClass(
  "gsList",
  slots = c(
    genesetsDataframe = "data.frame",
    enrichmentScore = "data.frame",
    #source = "ANY",
    cluster = "ANY",
    labels = "list"
  )
)
