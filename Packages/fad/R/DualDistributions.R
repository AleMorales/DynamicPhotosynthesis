#' @include DualVector.R
#' @include DualVectorMath.R
NULL



# Normal distribution ---------------------------------------------------------------------------------------------
setGeneric("dnorm")
dnorm_fun = function(x, mean, sd, log) {
  f = 1/sqrt(2*sd^2*pi)*exp(-(x - mean)^2/(2*sd^2))
  if(missing(log) | !log) {
    return(f)
  } else {
    return(log(f))
  }
}
#' @export
setMethod("dnorm", c("ANY", "DualVector", "DualVector", "logical"), dnorm_fun)
#' @export
setMethod("dnorm", c("ANY", "DualVector", "DualVector", "missing"), dnorm_fun)
#' @export
setMethod("dnorm", c("ANY", "DualVector", "ANY", "logical"), dnorm_fun)
#' @export
setMethod("dnorm", c("ANY", "ANY", "DualVector", "logical"), dnorm_fun)
#' @export
setMethod("dnorm", c("ANY", "ANY", "DualVector", "missing"), dnorm_fun)
#' @export
setMethod("dnorm", c("ANY", "DualVector", "ANY", "missing"), dnorm_fun)
