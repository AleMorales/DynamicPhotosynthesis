# For vector-valued functions


# Finite difference -----------------------------------------------------------------------------------------------
#' @export
finite_jac = function(ofun, method = "simple", h = 1e-4, side = 1, d = 0.0001,
                       zero.tol = sqrt(.Machine$double.eps/7e-7), r = 4, v = 2) {
  switch(method,

         simple = function(pars, ...) numDeriv::jacobian(ofun, pars, method = "simple", method.args = list(eps = h),
                                                     side = if(length(side) == 1) rep(side, length(pars)) else side, ...),

         Richardson = function(pars, ...) numDeriv::jacobian(ofun, pars, method = "Richardson", method.args = list(eps = h,
                                                         d = d, zero.tol = zero.tol, r = r, v = v), ...),

         central = function(pars, ...) {
           f = ofun(pars, ...)
           jac = matrix(nrow = length(f), ncol = length(pars))
           for(i in seq_along(pars)) {
             mod_pars = pars
             mod_pars[i] = mod_pars[i] - h
             fb = ofun(mod_pars, ...)
             mod_pars[i] = mod_pars[i] + 2*h
             ff = ofun(mod_pars, ...)
             jac[,i] = (ff - fb)/(2*h)
           }
           return(jac)
         },

         default = stop("Method ", method, " not supported.")
  )
}


# Complex step ----------------------------------------------------------------------------------------------------
#' @export
complex_jac = function(ofun) {
  function(pars, ...) {
    numDeriv::jacobian(ofun, pars, method = "complex", ...)
  }
}


# Symbolic --------------------------------------------------------------------------------------------------------
#' @export
symbolic_jac = function(ofun, parnames, ...) {
  f = Deriv::Deriv(ofun, parnames, ...)
  function(pars, ...) {
    out = f(pars, ...)
    matrix(unname(out), ncol = length(pars), byrow = TRUE)
  }
}


# Dual numbers ----------------------------------------------------------------------------------------------------
#' @export
dual_jac = function(ofun) {
  require(fad)
  function(pars, ...) fad::grad(ofun(fad::dualvector(pars, diag(length(pars))), ...))
}
