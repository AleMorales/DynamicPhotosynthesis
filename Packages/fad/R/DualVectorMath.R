#' @include DualVector.R
NULL


# Binary operators ------------------------------------------------------------------------------------------------

# Addition
#' @export
setMethod("+", c("DualVector", "DualVector"), function(e1, e2) {
  val = e1@val + e2@val
  grad = e1@grad + e2@grad
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("+", c("DualVector", "numeric"), function(e1, e2) {
  val = e1@val + e2
  grad = e1@grad
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("+", c("numeric", "DualVector"), function(e1, e2) {
  val = e1 + e2@val
  grad = e2@grad
  e2@val = val
  e2@grad = grad
  e2
})


# Substraction
#' @export
setMethod("-", c("DualVector"), function(e1) {
  val = -e1@val
  grad = -e1@grad
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("-", c("DualVector", "DualVector"), function(e1, e2) {
  val = e1@val - e2@val
  grad = e1@grad - e2@grad
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("-", c("DualVector", "numeric"), function(e1, e2) {
  val = e1@val - e2
  grad = e1@grad
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("-", c("numeric", "DualVector"), function(e1, e2) {
  val = e1 - e2@val
  grad = -e2@grad
  e2@val = val
  e2@grad = grad
  e2
})


# Multiplication
#' @export
setMethod("*", c("DualVector", "DualVector"), function(e1, e2) {
  val = e1@val*e2@val
  e1_grad = matrix(t(e1@grad), ncol = ncol(e1@grad), nrow = length(val), byrow = TRUE)
  e2_grad = matrix(t(e2@grad), ncol = ncol(e2@grad), nrow = length(val), byrow = TRUE)
  grad = e1_grad*e2@val + e1@val*e2_grad
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("*", c("DualVector", "numeric"), function(e1, e2) {
  val = e1@val*e2
  e1_grad = matrix(t(e1@grad), ncol = ncol(e1@grad), nrow = length(val), byrow = TRUE)
  grad = e1_grad*e2
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("*", c("numeric", "DualVector"), function(e1, e2) {
  val = e1*e2@val
  e2_grad = matrix(t(e2@grad), ncol = ncol(e2@grad), nrow = length(val), byrow = TRUE)
  grad =  e2_grad*e1
  e2@val = val
  e2@grad = grad
  e2
})

# Division
#' @export
setMethod("/", c("DualVector", "DualVector"), function(e1, e2) {
  val = e1@val/e2@val
  e1_grad = matrix(t(e1@grad), ncol = ncol(e1@grad), nrow = length(val), byrow = TRUE)
  e2_grad = matrix(t(e2@grad), ncol = ncol(e2@grad), nrow = length(val), byrow = TRUE)
  grad = (e1_grad*e2@val - e1@val*e2_grad)/e2@val^2
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("/", c("DualVector", "numeric"), function(e1, e2) {
  val = e1@val/e2
  e1_grad = matrix(t(e1@grad), ncol = ncol(e1@grad), nrow = length(val), byrow = TRUE)
  grad = e1_grad/e2
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("/", c("numeric", "DualVector"), function(e1, e2) {
  val = e1/e2@val
  e2_grad = matrix(t(e2@grad), ncol = ncol(e2@grad), nrow = length(val), byrow = TRUE)
  grad = -e1*e2_grad/e2@val^2
  e2@val = val
  e2@grad = grad
  e2
})


# Power
#' @export
setMethod("^", c("DualVector", "DualVector"), function(e1, e2) {
  val = e1@val^e2@val
  grad = e1@grad*e2@val*(e1@val^(e2@val - 1)) +  e2@grad*(e1@val^e2@val)*log(e1@val)
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("^", c("DualVector", "numeric"), function(e1, e2) {
  val = e1@val^e2
  grad = e1@grad*e2*(e1@val^(e2 - 1))
  e1@val = val
  e1@grad = grad
  e1
})
#' @export
setMethod("^", c("numeric", "DualVector"), function(e1, e2) {
  val = e1^e2@val
  grad = e2@grad*(e1^e2@val)*log(e1)
  e1@val = val
  e1@grad = grad
  e2
})

#' @export
setMethod("abs", c("DualVector"), function(x) {ifelse(x < 0, -x, x)})


# Power functions -------------------------------------------------------------------------------------------------

#' @export
setMethod("exp", c("DualVector"), function(x) {
  val = exp(x@val)
  grad = exp(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("expm1", c("DualVector"), function(x) {
  val = expm1(x@val)
  grad = expm1(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})

#' @export
setMethod("log", c("DualVector"), function(x) {
  val = log(x@val)
  grad = x@grad/x@val
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("log1p", c("DualVector"), function(x) {
  val = log1p(x@val)
  grad = x@grad/(1 + x@val)
  x@val = val
  x@grad = grad
  x
})

#' @export
setMethod("sqrt", c("DualVector"), function(x) {x^0.5})


# Trigonometry ----------------------------------------------------------------------------------------------------

#' @export
setMethod("sin", c("DualVector"), function(x) {
  val = sin(x@val)
  grad = cos(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("sinpi", c("DualVector"), function(x) {
  val = sinpi(x@val)
  grad = pi*cospi(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("cos", c("DualVector"), function(x) {
  val = cos(x@val)
  grad = -sin(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("cospi", c("DualVector"), function(x) {
  val = cospi(x@val)
  grad = -pi*sinpi(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("tan", c("DualVector"), function(x) {
  val = tan(x@val)
  grad = (tan(x@val)^2 + 1)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("tanpi", c("DualVector"), function(x) {
  val = tanpi(x@val)
  grad = pi/cospi(x@val)^2*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("asin", c("DualVector"), function(x) {
  val = asin(x@val)
  grad = 1/(sqrt(1 - x@val^2))*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("atan", c("DualVector"), function(x) {
  val = atan(x@val)
  grad = 1/(sqrt(1 + x@val^2))*x@grad
  x@val = val
  x@grad = grad
  x
})



# Hyperbolic trigonometry -----------------------------------------------------------------------------------------


#' @export
setMethod("sinh", c("DualVector"), function(x) {
  val = sinh(x@val)
  grad = cosh(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("cosh", c("DualVector"), function(x) {
  val = cosh(x@val)
  grad = sinh(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("tanh", c("DualVector"), function(x) {
  val = tanh(x@val)
  grad = (1 - tanh(x@val)^2)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("acosh", c("DualVector"), function(x) {
  val = acosh(x@val)
  grad = 1/sqrt(x@val^2 - 1)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("asinh", c("DualVector"), function(x) {
  val = asinh(x@val)
  grad = 1/sqrt(x@val^2 + 1)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("atanh", c("DualVector"), function(x) {
  val = atanh(x@val)
  grad = 1/(1 - x@val^2)*x@grad
  x@val = val
  x@grad = grad
  x
})
setGeneric("atan2")
#' @export
setMethod("atan2", c("DualVector", "DualVector"), function(y, x) {
  val = atan2(y@val, x@val)
  grad = x@val/(x@val^2 + y@val^2)*y@grad - y@val/(x@val^2 + y@val^2)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("atan2", c("DualVector", "numeric"), function(y, x) {
  val = atan2(y@val, x)
  grad = x/(x^2 + y@val^2)*y@grad
  y@val = val
  y@grad = grad
  y
})
#' @export
setMethod("atan2", c("numeric", "DualVector"), function(y, x) {
  val = atan2(y, x@val)
  grad = -y/(x@val^2 + y^2)*x@grad
  x@val = val
  x@grad = grad
  x
})

# Special functions -----------------------------------------------------------------------------------------------
#' @export
setMethod("gamma", c("DualVector"), function(x) {
  val = gamma(x@val)
  grad = gamma(x@val)*digamma(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("lgamma", c("DualVector"), function(x) {
  val = lgamma(x@val)
  grad = digamma(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("digamma", c("DualVector"), function(x) {
  val = digamma(x@val)
  grad = trigamma(x@val)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("trigamma", c("DualVector"), function(x) {
  val = trigamma(x@val)
  grad = psigamma(x@val, 2L)*x@grad
  x@val = val
  x@grad = grad
  x
})
setGeneric("psigamma")
#' @export
setMethod("psigamma", c("DualVector", "numeric"), function(x, deriv) {
  val = psigamma(x@val, deriv)
  grad = psigamma(x@val, 1L + deriv)*x@grad
  x@val = val
  x@grad = grad
  x
})
#' @export
setMethod("psigamma", c("DualVector", "integer"), function(x, deriv) {
  val = psigamma(x@val, deriv)
  grad = psigamma(x@val, 1L + deriv)*x@grad
  x@val = val
  x@grad = grad
  x
})



