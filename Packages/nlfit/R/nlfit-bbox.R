#' @include generics.R
NULL

# Fit model to data assuming Gaussian likelihood.
# model, error, jac, jac_error, ineq, ineq_grad, eq and eq_grad should all be closures
nlfit_bbox = function(
  model,
  obs,
  error = NULL,
  start,
  fixed = NULL,
  lower = NULL,
  upper = NULL,
  algorithm = "nelder",
  options = NULL,
  jac = "central",
  jac_error = "central",
  ineq = NULL,
  ineq_grad = NULL,
  eq = NULL,
  eq_grad = NULL,
  trace = FALSE,
  eval = FALSE) {

  # Add start and lower
  if(is.null(lower)) lower = rep(-Inf, length(start))
  if(is.null(upper)) upper = rep(Inf, length(start))

  # If jac is missing construct it using central differences
  if(is.null(jac))
    jacfun = diff::finite_jac(model, method = "central")
  else
    jacfun = jac

  # If error function is provided but its Jacobian is missing, use central differences as approximation
  merror = is.null(error)
  if(!merror) {
    if(is.character(jac_error)) {
      jac_error_fun = switch(jac_error,
                             central = diff::finite_jac(error, method = "central"),
                             simple = diff::finite_jac(error, method = "simple"),
                             complex = diff::complex_jac(error),
                             symbolic = {jactext = deparse(diff::symbolic_jac(error[[3]], names(start)), width.cutoff = 500)
                             jactext[max(length(jactext) - 1, 1)] = sub("c", replacement = "cbind", x = jactext[max(length(jactext) - 1,1)], fixed = TRUE)
                             eval(parse(text =  paste0("function() with(data, function(pars) with(as.list(pars),", paste0(jactext, collapse = "\n"), "))")))()},
                             stop("jac of type ", jac_error, " not available")
      )
    } else {
      jac_error_fun = jac_error
    }
  }

  # Construct optimization function from individual functions.
  # If merror, treat sd as nuisance parameter (= sd(residuals))
  NLL = function(pars) {
    mu = model(pars)
    if(merror) {
      n = length(mu)
      sigma = sqrt(sum((obs - mu)^2, na.rm = TRUE)/(n - 1))
    } else {
      sigma = error(pars)
    }
    -sum(-log(sqrt(2*pi)) - log(sigma) - (obs - mu)^2/(2*sigma^2), na.rm = TRUE)
  }

  # Construct gradient.
  # If merror, treat sd as nuisance parameter (= sd(residuals))
  NLL_grad = function(pars) {
    mu = model(pars)
    gmu = jacfun(pars)
    if(merror) {
      n = length(mu)
      s = sqrt(sum((obs - mu)^2, na.rm = TRUE)/(n - 1))
      gs = 0.5*sum(-2*(obs - mu)*gmu, na.rm = TRUE)/((n - 1)*s)
    } else {
      s = error(pars)
      gs = jac_error_fun(pars)
    }
    # Make sure each measurement point has an associated s
    if(length(s) == 1) {
      s = rep(s,times = length(mu))
      gs = matrix(rep(as.numeric(gs), times = length(mu)), nrow = length(mu), ncol = ncol(gs), byrow = TRUE)
    }
    colSums(1/s*gs - 1/s^3*(obs - mu)^2*gs - (obs - mu)/s^2*gmu, na.rm = TRUE)
  }


  # Minimize the ofun using mlefit
  result = mlefit(NLL = NLL, start = start, fixed = fixed, lower = lower, upper = upper, algorithm = algorithm,
                  options = options, grad = NLL_grad, ineq = ineq, ineq_grad = ineq_grad,
                  eq = eq, eq_grad = eq_grad, trace = trace, eval = eval)

  # Store the model as a callable function
  result$model = model

  # Additional inputs required to rerun fit (or model evaluation)
  result$lower = lower
  result$upper = upper
  result$options$parscale = options$parscale

  # Store original obs and fitted (take into account fixed parameters)
  result$obs = obs
  if(!is.null(fixed)) {
    all_par = start
    imod = seq_along(start)[-which(names(start) %in% fixed)]
    all_par[imod] = result$par
    result$fitted = model(all_par)
    result$all_par = all_par
  } else {
    result$fitted = model(result$par)
    result$all_par = result$par
  }

  # Store the all the information to refit the model again (for profiling)
  result$error = error
  result$jac = jac
  result$jac_error = jac_error
  result$fixed = fixed


  # Assign nlfit class
  class(result) = c("nlfit_bbox", "nlfit")

  result

}


# Make predictions
#' @export
#' @method predict nlfit_bbox
predict.nlfit_bbox = function(object) {
  object$model(object$par)
}


# Refit the model with possibility of changing any settings or inputs
#' @export
#' @method update nlfit_bbox
update.nlfit_bbox = function(object, model = object$model, error = object$error, obs = object$obs,
                             start = object$all_par, fixed = object$fixed, lower = object$lower, upper = object$upper,
                             algorithm = object$algorithm, jac = object$jac, jac_error = object$jac_error,
                             ineq = object$ineq, ineq_grad = object$ineq_grad, eq = object$eq, eq_grad = object$eq_grad,
                             options, trace = FALSE, eval = FALSE) {

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
  nlfit_bbox(model = model, error = error, obs = obs, start = start,
             fixed = fixed, lower = valid_lower, upper = valid_upper, algorithm =  algorithm,
             options = object$options, jac =  jac, jac_error =  jac_error,
             ineq =  ineq, ineq_grad =  ineq_grad, eq =  eq, eq_grad =  eq_grad,
             trace = trace, eval = eval)
}
