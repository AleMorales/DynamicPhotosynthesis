#' @include generics.R profiler.R hessian.R slicer.R
NULL


#' @export
mlefit = function(NLL, start, fixed = NULL,
                  lower = rep(-Inf, length(start)),
                  upper = rep(Inf, length(start)),
                  algorithm = "nelder",
                  options = NULL,
                  grad = NULL,
                  ineq = NULL,
                  ineq_grad = NULL,
                  eq = NULL,
                  eq_grad = NULL,
                  trace = FALSE,
                  eval = FALSE) {
  # Check inputs
  if (!is.function(NLL))
    stop("The argument NLL must be a function.")
  if(is.null(names(upper)) && (length(start) != length(upper)))
    stop("The length of start and upper must be the same (unless using named parameters in upper).")
  if(is.null(names(lower)) && (length(start) != length(lower)))
    stop("The length of start and lower must be the same (unless using named parameters in lower).")
  if (!is.character(algorithm))
    stop("The argument algorithm must be a string.")
  if (!is.null(options) && !is.list(options))
    stop("The argument options must be a list")

  # Parameters within start can be passed unchanged to the model as specified in argument fixed
  # fixed can be a list of parameter names or positions within argument start
  if(is.null(fixed)) {
    imod = seq_along(start)
  } else {
    if(is.character(fixed)) {
      if(!all(fixed %in% names(start))) stop("Argument fixed must refer to named parameters in argument start")
      imod = seq_along(start)[-which(names(start) %in% fixed)]
    } else if(is.numeric(fixed)) {
      if(any(fixed > length(start) || fixed < 1)) stop("Argument fixed must refer to positions within argument start")
      imod = seq_along(start)[-fixed]
    } else {
      stop("Argument fixed must be a numeric or character vector")
    }
  }

  # Create valid lower and upper arguments
  # They can be unnamed (they already valid) or named (need to create new lower/upper with length of start)
  if(is.null(names(lower))) {
    valid_lower = lower
  } else {
    valid_lower = rep(-Inf, length(start))
    names(valid_lower) = names(start)
    valid_lower[names(lower)] = unname(lower)
  }
  if(is.null(names(upper))) {
    valid_upper = upper
  } else {
    valid_upper = rep(Inf, length(start))
    names(valid_upper) = names(start)
    valid_upper[names(upper)] = unname(upper)
  }


  # Create wrapper around all functions to ensure that fixed parameters are passed correctly
  ofun = function(par) {.x = start; .x[imod] = par; NLL(.x)}
  grad_fun = function(par) {.x = start; .x[imod] = par; grad(.x)[imod]}
  ineq_fun = if(!is.null(ineq)) function(par) {.x = start; .x[imod] = par;  ineq(.x)} else NULL
  eq_fun = if(!is.null(eq)) function(par) {.x = start; .x[imod] = par; eq(.x)} else NULL
  ineq_grad_fun = if(!is.null(ineq_grad)) function(par) {.x = start; .x[imod] = par; ineq_grad(.x)} else NULL
  eq_grad_fun = if(!is.null(eq_grad)) function(par) {.x = start; .x[imod] = par; eq_grad(.x)} else NULL

  # Modified parscale
  if(!is.null(options$parscale))
    options$parscale = options$parscale[imod]

  # Run optimization if eval = FALSE
  if(!eval) {
    result = RcppNLopt::nlopt(x0 = start[imod], ofun = ofun, lb = valid_lower[imod], ub = valid_upper[imod],
                              algorithm = algorithm,
                              options = options, grad = grad_fun, ineq = ineq_fun, ineq_grad = ineq_grad_fun,
                              eq = eq_fun, eq_grad = eq_grad_fun, trace = trace)
  } else {
    result = list(min = ofun(start[imod]), par = start[imod], neval = 1, algorithm = NA, message = NA)
  }
  # Adapt names to the mlefit terminology
  names(result)[1] = "NLL"

  # Store inputs to construct modelling problem again (because things like imod may change)
  result$NLL_fun = ofun
  result$grad = grad
  result$ineq = ineq
  result$ineq_grad = ineq_grad
  result$eq = eq
  result$eq_grad = eq_grad

  # Store lower and upper boundaries
  result$lower = lower
  result$upper = upper

  # Assign names to parameters
  names(result$par) = names(start[imod])
  result$all_par = start
  result$all_par[imod] = result$par

  # Store parameters that were fixed
  result$fixed = fixed

  # Store options used int the simulations
  result$options = options

  # mlefit class to overload methods
  class(result) = "mlefit"

  # Return
  result
}


# Basic mlefit methods --------------------------------------------------------------------------------------------

# Extract coefficients
#' @export
#' @method coef mlefit
coef.mlefit = function(object, ...) {
  object$par
}


# Default print of the object
#' @export
#' @method print mlefit
print.mlefit = function(x, ...) {
  cat("Coefficients:\n")
  print(coef(x))
  cat("\nNLL: ")
  cat(x$NLL, "\n")
}


# Print diagnostics of the optimization
#' @export
#' @method diagnostics mlefit
diagnostics.mlefit = function(x) {
  print(x[c("algorithm", "status", "message", "nevals")])
}


# Refit the model with possibility of changing any settings or inputs
#' @export
#' @method update mlefit
update.mlefit = function(object, NLL = object$NLL_fun, start = object$all_par,
                         fixed = object$fixed, lower = object$lower, upper = object$upper,
                         algorithm = object$algorithm, grad = object$grad, ineq = object$ineq,
                         ineq_grad = object$ineq_grad, eq = object$eq, eq_grad = object$eq_grad,
                          options, trace = FALSE, eval = FALSE, ...) {

  # Allows to specify partial changes in options
  if(!missing(options))
  for(i in names(options)) object$options[[i]] = options[[i]]

  # Allows to specify partial changes to lower and upper (include parameters that were not bound before)
  if(is.null(names(lower))) {
    if(length(lower) != length(object$lower)) stop("Names missing from lower and length is not valid")
    valid_lower = lower
  } else {
    if(!all(names(lower) %in% names(object$all_par))) stop("Name within lower does not correspond to a parameter of the model.")
    valid_lower = object$lower
    if(is.null(names(valid_lower))) names(valid_lower) = names(object$all_par)
    valid_lower[names(lower)] = unname(lower)
  }
  if(is.null(names(upper))) {
    if(length(upper) != length(object$upper)) stop("Names missing from upper and length is not valid")
    valid_upper = upper
  } else {
    if(!all(names(upper) %in% names(object$all_par))) stop("Name within upper does not correspond to a parameter of the model.")
    valid_upper = object$upper
    if(is.null(names(valid_upper))) names(valid_upper) = names(object$all_par)
    valid_upper[names(upper)] = unname(upper)
  }

  # Refit the model
  mlefit(NLL = NLL, start = start, fixed = fixed, lower = valid_lower,
         upper = valid_upper, algorithm = algorithm, options = object$options,
         grad = grad, ineq = ineq, ineq_grad = ineq_grad, eq = eq,
         eq_grad = eq_grad, trace = trace, eval = eval, ...)
}

# Calculate the likelihood profiles for different parameters
#' @export
#' @method profile mlefit
profile.mlefit = profiler

# Calculate the likelihood profiles for different parameters
#' @export
#' @method slice mlefit
slice.mlefit = slicer

# Calculate confidence interval with different possible methods
#' @export
#' @method confint mlefit
confint.mlefit = function(object, type = "profile", parm = names(object$par), level = 0.95, ...) {
  switch(type,
    profile = {profiles = profile(object, parm = parm, level = level, ...); confint(profiles)},
    slice = {slices = slice(object, parm = parm, level = level, ...); confint(slices)},
    hessian = calc_confint_hessian(object, parm, level, ...),
    bootstrap = confint_bootstrap(object, parm, level, ...),
    stop("Method ", type, " for calculating confidence intervals not supported.")
  )
}


# Calculate the variance covariance matrix of the model
#' @export
#' @method vcov mlefit
vcov.mlefit = calc_vcov
