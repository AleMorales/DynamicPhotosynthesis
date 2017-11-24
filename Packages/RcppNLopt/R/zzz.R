#' @useDynLib RcppNLopt
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("RcppNLopt", libpath)
}
