#' @include DualVector.R
NULL


# Reduction operations --------------------------------------------------------------------------------------------

#' @export
setMethod("sum", c("DualVector"), function(x) {
  x@val = sum(x@val)
  x@grad = matrix(colSums(x@grad), ncol = ncol(x@grad))
  x
})
#' @export
setMethod("mean", c("DualVector"), function(x) {
  x@val = mean(x@val)
  x@grad = matrix(colSums(x@grad)/length(x@val), ncol = ncol(x@grad))
  x
})
#' @export
setMethod("median", c("DualVector"), function(x) {
  n = length(x)
  if(n %% 2) {
    x@val = median(x@val)
    x@grad = matrix(x@grad[n %/% 2 + 1, ], ncol = ncol(x@grad))
  } else {
    x@val = median(x@val)
    x@grad = matrix(colMeans(x@grad[c(n %/% 2 + 1, n %/% 2), , drop = FALSE]), ncol = ncol(x@grad))
  }
  x
})


#' @export
setMethod("max", c("DualVector"), function(x) {
  i = which.max(x@val)
  x@val = x@val[i]
  x@grad = x@grad[i,]
  x
})
#' @export
setMethod("min", c("DualVector"), function(x) {
  i = which.min(x@val)
  x@val = x@val[i]
  x@grad = x@grad[i,]
  x
})







#' @export
setMethod("cumsum", c("DualVector"), function(x) {
  x@val = cumsum(x@val)
  x@grad = matrix(cumsum(x@grad), ncol = ncol(x@grad))
  x
})

#' @export
setMethod("cumprod", c("DualVector"), function(x) {
  stop("Method not implemented")
})

#' @export
setMethod("cummax", c("DualVector"), function(x) {
  stop("Method not implemented")
})

#' @export
setMethod("cummin", c("DualVector"), function(x) {
  stop("Method not implemented")
})

