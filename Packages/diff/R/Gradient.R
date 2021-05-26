# For multivariate and scalar functions
# Multiple inputs must be passed through the first argument


# Finite difference -----------------------------------------------------------------------------------------------
#' @export
finite_grad = function(ofun, method = "simple", h = 1e-4, side = +1, d = 0.0001,
                       zero.tol = sqrt(.Machine$double.eps/7e-7), r = 4, v = 2) {
  switch(method,

         simple = function(pars, ...) numDeriv::grad(ofun, pars, method = "simple", method.args = list(eps = h),
                                                     side = if(length(side) == 1) rep(side, length(pars)) else side, ...),

         Richardson = function(pars, ...) numDeriv::grad(ofun, pars, method = "Richardson", method.args = list(eps = h,
                                                         d = d, zero.tol = zero.tol, r = r, v = v), ...),

         central = function(pars, ...) {
           gr = numeric(length(pars))
           for(i in seq_along(pars)) {
             mod_pars = pars
             mod_pars[i] = mod_pars[i] - h
             fb = ofun(mod_pars, ...)
             mod_pars[i] = mod_pars[i] + 2*h
             ff = ofun(mod_pars, ...)
             gr[i] = (ff - fb)/(2*h)
           }
           return(gr)
         },

         default = stop("Method ", method, " not supported.")
  )
}


# Complex step ----------------------------------------------------------------------------------------------------
#' @export
complex_grad = function(ofun) {
  function(pars, ...) {
    numDeriv::grad(ofun, pars, method = "complex", ...)
  }
}


# Symbolic --------------------------------------------------------------------------------------------------------
#' @export
symbolic_grad = function(ofun, parnames, ...) {
   Deriv::Deriv(ofun, parnames, ...)
}


# Dual numbers ----------------------------------------------------------------------------------------------------
#' @export
dual_grad = function(ofun) {
  require(fad)
  function(pars, ...) as.numeric(fad::grad(ofun(fad::dualvector(pars, diag(length(pars))), ...)))
}

