#' @include generics.R
NULL

# Fit model to data assuming Gaussian likelihood.
# Model, error, ineq and eq specified as formulae
# Jacobians and gradients specified as functions
nlfit_formula = function(
  model,
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
  data,
  trace = FALSE,
  eval = FALSE) {

  # Add start and lower
  if(is.null(lower)) lower = rep(-Inf, length(start))
  if(is.null(upper)) upper = rep(Inf, length(start))

  # Convert formula into model
  model_formula = model

  # Extract vector of observations
  if(inherits(data, "tbl"))
    obs = data[,deparse(model[[2]])][[1]]
  else
    obs = data[,deparse(model[[2]])]

  model = eval(parse(text = paste0("function() with(data, function(pars) with(as.list(pars),", deparse(model[[3]], width.cutoff = 500), "))")))()

  # Generate Jacobian function of the model if missing
  if(is.character(jac)) {
    jacfun = switch(jac,
                   central = diff::finite_jac(model, method = "central"),
                   simple = diff::finite_jac(model, method = "simple"),
                   complex = diff::complex_jac(model),
                   symbolic = {jactext = deparse(diff::symbolic_jac(model_formula[[3]], names(start)), width.cutoff = 500)
                               jactext[max(length(jactext) - 1, 1)] = sub("c", replacement = "cbind", x = jactext[max(length(jactext) - 1,1)], fixed = TRUE)
                               eval(parse(text =  paste0("function() with(data, function(pars) with(as.list(pars),", paste0(jactext, collapse = "\n"), "))")))()},
                   stop("jac of type ", jac, " not available")
    )
  } else {
    jacfun = jac
  }

  # If error formula is provided construct function and jacobian
  merror = is.null(error)
  if(!merror) {
    errorfun = eval(parse(text = paste0("function() with(data, function(pars) with(as.list(pars),", deparse(error[[3]], width.cutoff = 500), "))")))()
    if(is.character(jac_error)) {
      jac_error_fun = switch(jac_error,
                             central = diff::finite_jac(errorfun, method = "central"),
                             simple = diff::finite_jac(errorfun, method = "simple"),
                             complex = diff::complex_jac(errorfun),
                             symbolic = {jactext = deparse(diff::symbolic_jac(error[[3]], names(start)), width.cutoff = 500)
                             jactext[max(length(jactext) - 1, 1)] = sub("c", replacement = "cbind", x = jactext[max(length(jactext) - 1,1)], fixed = TRUE)
                             eval(parse(text =  paste0("function() with(data, function(pars) with(as.list(pars),", paste0(jactext, collapse = "\n"), "))")))()},
                             stop("jac of type ", jac_error, " not available")
      )
    } else {
      jac_error_fun = jac_error
    }
  }


  # Construct optimization function from individual functions, depending on merror
  NLL = function(pars) {
    mu = model(pars)
    if(merror) {
      n = length(mu)
      sigma = sqrt(sum((obs - mu)^2, na.rm = TRUE)/(n - 1))
    } else {
      sigma = errorfun(pars)
    }
    -sum(-log(sqrt(2*pi)) - log(sigma) - (obs - mu)^2/(2*sigma^2), na.rm = TRUE)
  }

  # Construct gradient of NLL
  NLL_grad = function(pars) {
    mu = model(pars)
    gmu = jacfun(pars)
    if(merror) {
      n = length(mu)
      s = sqrt(sum((obs - mu)^2, na.rm = TRUE)/(n - 1))
      gs = 0.5*sum(-2*(obs - mu)*gmu, na.rm = TRUE)/((n - 1)*s)
    } else {
      s = errorfun(pars)
      gs = jac_error_fun(pars)
    }
    if(length(s) == 1) {
      s = rep(s,times = length(mu))
      gs = matrix(rep(as.numeric(gs), times = length(mu)), nrow = length(mu), ncol = ncol(gs), byrow = TRUE)
    }
    colSums(1/s*gs - 1/s^3*(obs - mu)^2*gs - (obs - mu)/s^2*gmu, na.rm = TRUE)
  }

  # Create function for ineq
  if(!is.null(ineq))
    ineq_fun = eval(parse(text = paste0("function() with(data, function(pars) with(as.list(pars),", deparse(ineq[[2]], width.cutoff = 500), "))")))()
  else
    ineq_fun = ineq

  # Minimize the ofun using mlefit
  result = mlefit(NLL = NLL, start = start, fixed = fixed, lower = lower, upper = upper, algorithm = algorithm,
                  options = options, grad = NLL_grad, ineq = ineq_fun, ineq_grad = ineq_grad,
                  eq = eq, eq_grad = eq_grad, trace = trace, eval = eval)

  # Store the model_formula (required to perform predictions)
  result$model_formula = model_formula

  # Additional inputs required to rerun fit (or model evaluation)
  result$lower = lower
  result$upper = upper
  result$options$parscale = options$parscale

  # Store the data
  result$data = data

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

  #

  # Assign nlfit class
  class(result) = c("nlfit_formula", "nlfit")

  result

}

# Make predictions
#' @export
#' @method predict nlfit_formula
predict.nlfit_formula = function(object, newdata = object$data) {
  model = object$model_formula
  model = eval(parse(text = paste0("function() with(newdata, function(pars) with(as.list(pars),", deparse(model[[3]], width.cutoff = 500), "))")))()
  model(object$par)
}

# Refit the model with possibility of changing any settings or inputs
#' @export
#' @method update nlfit_formula
update.nlfit_formula = function(object, model = object$model, error = object$error, data = object$data,
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
  nlfit_formula(model = model, error = error, start = start, fixed = fixed, lower = valid_lower,
         upper = valid_upper, algorithm = algorithm, options = object$options,
         jac = jac, jac_error = jac_error, data = data, ineq = ineq, ineq_grad = ineq_grad,
         eq = eq, eq_grad = eq_grad, trace = trace, eval = eval)
}
