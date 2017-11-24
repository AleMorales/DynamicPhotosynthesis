calc_vcov = function (object, parm = names(object$par), hessian = "Richardson", ...) {

  if(is.character(hessian)) {

    # Wrapper to evaluate NLL for different parameter values
    f = function(par, ...) {
      start = object$all_par
      start[parm] = par
      update(object, start = start, eval = TRUE, ...)$NLL
    }

    # Calculate hessian matrix for selected parameters
    hess = switch(hessian,
                  Richardson = diff::finite_hess(f, ...),
                  complex = diff::complex_hess(f),
                  stop("The hessian method ", hessian, " is not available."))
    hessian_mle = hess(object$par[parm], ...)

  } else if(is.function(hessian)) {
    hessian_mle = hessian(object$par[parm], ...)
  } else {
    stop("Argument hessian must be a string or a function")
  }

  # Calculate vcov matrix (we do not take negative as we used Negative Log-Likelihood)
  solve(hessian_mle)
}

# Calculate confidence intervals by calculating the variance-covariance matrix
calc_confint_hessian = function(object, parm, level, hessian = "Richardson", ...) {

  # Variance-covariance matrix
  vcov = vcov(object, parm, hessian, ...)

  # Calculate standard error of each parameter
  se = sqrt(diag(vcov))

  # Coefficient between confidence interval at a given level and sigma
  upper = object$par[parm] + qnorm((1 - level)/2, 0, 1, lower.tail = FALSE)*se
  lower = object$par[parm] - qnorm((1 - level)/2, 0, 1, lower.tail = FALSE)*se

  # Return matrix with mle and confidence intervals
  out = data.frame(value = object$par[parm], lower = lower, upper = upper)
  row.names(out) = parm
  out

}
