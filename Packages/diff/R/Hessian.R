# For multivariate- functions


# Finite difference -----------------------------------------------------------------------------------------------
#' @export
finite_hess = function(ofun, h = 1e-4, d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-7), r = 4, v = 2, ...) {
  function(pars, ...) numDeriv::hessian(func = ofun, x = pars, method = "Richardson",
                                         method.args = list(eps = h, d = d, zero.tol = zero.tol, r = r, v = v), ...)
}

# Complex step ----------------------------------------------------------------------------------------------------
#' @export
complex_hess = function(ofun) {
  function(pars, ...) {
    numDeriv::hessian(func = ofun, x = pars, method = "complex", ...)
  }
}

# Symbolic --------------------------------------------------------------------------------------------------------
#' @export
symbolic_hess = function(ofun, parnames, ...) {
  f = Deriv::Deriv(ofun, parnames, nderiv = 2, ...)
  function(pars, ...) {
    out = f(pars, ...)
    matrix(unname(out), ncol = length(pars), byrow = TRUE)
  }
}
