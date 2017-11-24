#' @include nlopt_options.R
NULL

#' @export
nlopt = function(x0,
                 ofun,
                 lb = rep(-Inf, length(x0)),
                 ub = rep(Inf, length(x0)),
                 algorithm = "nelder",
                 options = NULL,
                 grad = NULL,
                 ineq = NULL,
                 ineq_grad = NULL,
                 eq = NULL,
                 eq_grad = NULL,
                 trace = FALSE,
                 ...) {

  # Check lengths of ub, lb and x0
  if(length(ub) != length(x0) || length(lb) != length(x0))
    stop("Lengths of x0, lb and ub must coincide.")

  # Construct valid options object with proper defaults
  nlopt_options = default_options(algorithm)
  if(!is.null(options)) {
    for(i in names(options)) {
      if(!(i %in% names(nlopt_options)))
        warning("Option ", i, " ignored.")
      nlopt_options[[i]] = options[[i]]
    }
    nlopt_options = update_nlopt_algorithm(nlopt_options)
  }

  # Check parscale has the right length
  if(length(nlopt_options$parscale) != 1 && length(nlopt_options$parscale) != length(x0))
    stop("Length of parscale and x0 must coincide.")

  # Adjust dimensions of the problem with parscale
  x0 = x0/nlopt_options$parscale
  ub = ub/nlopt_options$parscale
  lb = lb/nlopt_options$parscale

  # Closure for ofun
  names_args = names(list(...))
  #checkClosure(ofun, names_args, deparse(substitute(ofun)))
  names_x = names(x0)
  if(trace == "plot") {
    par(mfrow = c(1,1))
    cumout = rep(Inf, length = nlopt_options$maxiter)
    ITER = 0
    nlopt_ofun = function(x) {
      names(x) = names_x
      out = nlopt_options$fnscale*ofun(nlopt_options$parscale*x, ...)
      ITER <<- ITER + 1
      cumout[ITER] <<- min(cumout,out, na.rm = TRUE)
      init = max(1,ITER - 50)
      plot(init:ITER, cumout[init:ITER], t = "l")
      out
    }
  } else if(trace == FALSE) {
    nlopt_ofun = function(x) {
      names(x) = names_x
      nlopt_options$fnscale*ofun(nlopt_options$parscale*x, ...)
    }
  } else if(trace == TRUE) {
    nlopt_ofun = function(x) {
      names(x) = names_x
      out = nlopt_options$fnscale*ofun(nlopt_options$parscale*x, ...)
      cat("Parameters:\n")
      cat(nlopt_options$parscale*x, "\n")
      cat("Objective:\n")
      cat(out, "\n")
      out
    }
  }

  test = try(nlopt_ofun(x0))
  if(inherits(test, "try-error")) stop("Error when calling objective function")

  # Closures for constraints
  if(!is.null(ineq)) {
    #checkClosure(ineq, names_args, deparse(substitute(ineq)))
    nlopt_ineq = function(x) as.matrix(ineq(nlopt_options$parscale*x, ...))
    test = try(nlopt_ineq(x0))
    if(inherits(test, "try-error")) stop("Error when calling inequality constraint")
  } else {
    nlopt_ineq = NULL
  }
  if(!is.null(eq)) {
    #checkClosure(eq, names_args, deparse(substitute(eq)))
    nlopt_eq = function(x) as.matrix(eq(nlopt_options$parscale*x, ...))
    test = try(nlopt_eq(x0))
    if(inherits(test, "try-error")) stop("Error when calling equality constraint")
  } else {
    nlopt_eq = NULL
  }

  # Add closures for gradients when needed (i.e. non-df algorithm)
  if(!(algorithm %in% df_algorithms)) {
    # Gradient of ofun
    if(!is.null(grad)) {
      #checkClosure(grad, names_args, deparse(substitute(grad)))
      names_x = names(x0)
      nlopt_grad = function(x) {
        names(x) = names_x
        nlopt_options$fnscale*grad(nlopt_options$parscale*x, ...)
      }
      test = try(nlopt_grad(x0))
      if(inherits(test, "try-error")) stop("Error when calling gradient of objective function")
    } else {
      stop("Algorithm ", algorithm, " requires the gradient of the objective function (grad).")
    }
    # Gradient of ineq if present
    if(!is.null(ineq)) {
      if(!is.null(ineq_grad)) {
        #checkClosure(ineq_grad, names_args, deparse(substitute(ineq_grad)))
        names_x = names(x0)
        nlopt_ineq_grad = function(x) {
          names(x) = names_x
          as.matrix(ineq_grad(nlopt_options$parscale*x, ...))
        }
        test = try(nlopt_ineq_grad(x0))
        if(inherits(test, "try-error")) stop("Error when calling gradient of inequality")
      } else {
        stop("Algorithm ", algorithm, " requires the gradient of the inequality constraint function (ineq_grad).")
      }
    } else {
      nlopt_ineq_grad = NULL
    }
    # Gradient of eq if present
    if(!is.null(eq)) {
      if(!is.null(eq_grad)) {
        #checkClosure(eq_grad, names_args, deparse(substitute(eq_grad)))
        names_x = names(x0)
        nlopt_eq_grad = function(x) {
          names(x) = names_x
          as.matrix(eq_grad(nlopt_options$parscale*x, ...))
        }
        test = try(nlopt_eq_grad(x0))
        if(inherits(test, "try-error")) stop("Error when calling gradient of equality")
      }
      else
        stop("Algorithm ", algorithm, " requires the gradient of the equality constraint function (eq_grad).")
    } else {
      nlopt_eq_grad = NULL
    }
  } else {
    nlopt_grad = NULL
    nlopt_ineq_grad = NULL
    nlopt_eq_grad = NULL
  }

  # When required, check that tolerances are provided for all constraints
  nlopt_options = check_constraint("ineq", nlopt_ineq, x0, nlopt_options)
  nlopt_options = check_constraint("eq", nlopt_eq, x0, nlopt_options)

  # If constraints are present but not supported, need to use auglag
  nlopt_options = check_auglag(nlopt_options, nlopt_ineq, nlopt_eq)

  # Call Rcpp interface to nlopt
  result = .Call('RcppNLopt_wrap_nlopt', PACKAGE = 'RcppNLopt',
                 x0, lb, ub, nlopt_options,
                 nlopt_ofun, nlopt_grad, nlopt_ineq, nlopt_ineq_grad,
                 nlopt_eq, nlopt_eq_grad)

  # Add additional information to the result
  result$algorithm = algorithm
  #result$options = nlopt_options

  # Restore scaling
  result$par = result$par*nlopt_options$parscale
  result$min = result$min/nlopt_options$fnscale

  # S3 class for overloading methods
  class(result) = "nlopt"
  return(result)

}




# Check closures:
# fun = the function provided by the user
# names_args = names of the variables passed through "..."
# fname = name of the function passed to nlopt
checkClosure = function(fun, names_args, fname) {

  if(!is.function(fun)) stop(fname, " must be a function")
  # Get names of formal arguments and compare to names in "..."
  names_fargs = names(formals(fun)[-1])
  nfargs = length(names_fargs)
  nargs = length(names_args)
  # Check all possible mismatches
  if(nfargs == 0) {
    if(nargs > 0) stop(fname, " does not require extract inputs but you passed ", nargs, " extra inputs through nlopt")
  } else {
    if(nargs == 0) stop(fname, " requires extra arguments but none were passed.")
    # Gather arguments that are missing and not required
    missing = names_fargs[!(names_fargs %in% names_args)]
    notrequired = names_args[!(names_args %in% names_fargs)]
    if(length(missing) > 0) stop("The following arguments are required by ", fname, " but were NOT provided: ",
                                 paste0(missing, collapse = ", "))
    if(length(notrequired) > 0) stop("The following arguments are NOT required by ", fname, " but were provided: ",
                                     paste0(missing, collapse = ", "))
  }


}



# Check that length of eqatol and ineqatol coincides with length of contraints
check_constraint = function(type, fun, x0, nlopt_options) {
  if(!is.null(fun)) {
    test = fun(x = x0)
    m = length(test)
    # When using local_opts
    if(!is.null(nlopt_options$local_opts)) {
      m_user = length(nlopt_options$local_opts[[paste0(type,"atol")]])
      if(m != m_user)
        if(m_user == 1)
          nlopt_options$local_opts[[paste0(type,"atol")]] =
            rep(nlopt_options$local_opts[[paste0(type,"atol")]], m)
        else
          stop("Number of catol does not coincide with number of constraints")
        # When not using local_opts
    } else {
      m_user = length(nlopt_options[[paste0(type,"atol")]])
      if(m != m_user)
        if(m_user == 1)
          nlopt_options[[paste0(type,"atol")]] =
            rep(nlopt_options[[paste0(type,"atol")]], m)
        else
          stop("Number of catol does not coincide with number of constraints")
    }
  }
  # original or updated options
  return(nlopt_options)
}


# Check if auglag is needed and return an updated nlopt_options object
check_auglag = function(nlopt_options, nlopt_ineq, nlopt_eq) {

  # Not needed if both constraints are misisng
  if(is.null(nlopt_ineq) && is.null(nlopt_eq)) return(nlopt_options)

  # Some algorithm support constraints "natively"
  alg = nlopt_options$nlopt_algorithm
  if(alg %in% c("LD_SLSQP", "LN_COBYLA", "GN_ISRES")) return(nlopt_options)
  if(is.null(nlopt_eq) && alg == "LD_MMA") return(nlopt_options)

  # For the rest of cases, auglag is needed.
  if(nlopt_options$algorithm %in% df_algorithms)
    auglag = "LN_AUGLAG"
  else
    auglag = "LD_AUGLAG"
  nlopt_options$local_opts = nlopt_options[c("xatol", "xrtol", "fatol", "frtol", "maxiter", "maxtime")]
  nlopt_options$nlopt_algorithm = auglag
  nlopt_options
}
