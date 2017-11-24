# These are generics that are common to all classes in the package

#' Print technical information on the fitting of a model
#'
#' @param x Model object belonging to the classes \code{\link{mlefit}}, \code{bayfit} or \code{nlm}.
#'
#' @export
diagnostics = function(x) {
  UseMethod("diagnostics")
}


#' Calculate likelihood slice
#' @export
slice = function(object, ...) {
  UseMethod("slice")
}

