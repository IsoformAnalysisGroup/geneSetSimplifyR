setMethod("show", "gsList", function(object) {
  if( length( object@labels ) ) {
    print(
      stringr::str_c(
        'A gsList object with ',
        dplyr::n_distinct(object@enrichmentScore$pathway),
        ' gene sets from ',
        dplyr::n_distinct(object@enrichmentScore$source),
        ' source(s)',
        ' clustered into ',
        length(levels( object@cluster@active.ident )),
        ' labeled clusters.'
      )
    )
  } else if( !is.null(object@cluster) ) {
    print(
      stringr::str_c(
        'A gsList object with ',
        dplyr::n_distinct(object@enrichmentScore$pathway),
        ' gene sets from ',
        dplyr::n_distinct(object@enrichmentScore$source),
        ' source(s)',
        ' clustered into ',
        length(levels( object@cluster@active.ident )),
        ' clusters.'
      )
    )
  } else {
    print(
      stringr::str_c(
        'A gsList object with ',
        dplyr::n_distinct(object@enrichmentScore$pathway),
        ' gene sets from ',
        dplyr::n_distinct(object@enrichmentScore$source),
        ' source(s)'
      )
    )
  }

})
