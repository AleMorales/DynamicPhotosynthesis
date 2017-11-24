#' @importFrom mc2 mcmc
#' @importFrom coda mcmc.list
#' @include generics.R
NULL


#' Bayesian integration of a posterior
#'
#' @param lp Function that evaluates log-posterior given a vector of parameter values
#' @param start Initial values for the parameter values
#' @param algorithm A string indicating the MCMC algorithm to be used
#' @param options List of settings of the MCMC algorithm. See Details below.
#' @param cl Parallel cluster for running multiple chains in parallel
#'
#' @export
bayfit = function(lp, start, algorithm = "metropolis", options = list(), cl = NULL) {

  # Run mcmc to obtain sample of posterior distribution
  result = mcmc(lp, start, algorithm, options, cl)

  # Create the bayfit data structure
  if(is.null(dim(start))) {
    result = structure(.Data = list(sample = coda::mcmc(result, thin = attr(result, "thin")),
                                    lp = lp, start = start, algorithm = algorithm,
                                    options = options), class = "bayfit")
  } else {
    result = structure(.Data = list(sample = mcmc.list(lapply(result, function(x) coda::mcmc(x, thin = attr(x, "thin")))),
                                    lp = lp, start = start, algorithm = algorithm,
                                    options = options), class = "bayfit")
  }

  return(result)
}


#' @export
update.bayfit = function(object, start, algorithm = object$algorithm, options, cl = object$cl) {

  # Allows to specify partial changes in options
  if(!missing(options))
  for(i in names(options)) object$options[[i]] = options[[i]]

  # Start at a user-defined location or obtain it from stored sample
  if(missing(start)) {
    if(inherits(object$sample, "mcmc.list")) {
      nchains = length(object$sample)
      niter = nrow(object$sample[[1]])
      start = matrix(NA, ncol = ncol(object$sample[[1]]), nrow = nchains)
      for(i in 1:nchains) {
        start[i,] = object$sample[[i]][niter,]
      }
    } else {
      nchains = 1
      start = object$sample[nrow(object$sample),]
    }
  }

  # Generate new sample from posterior
  result = bayfit(object$lp, start, algorithm, object$options, cl)

  # Combine samples with the previous run
  if(nchains == 1) {
    result$sample = coda::mcmc(data = rbind(object$sample, result$sample), thin = attr(result$sample, "thin"))
    attr(result$sample, "thin") = attr(result$sample, "thin")
    attr(result$sample, "accept") = attr(result$sample, "accept")
  } else {
    for(i in 1:nchains) {
      result$sample[[i]] = coda::mcmc(rbind(object$sample[[i]], result$sample[[i]]), thin = attr(result$sample[[i]], "thin"))
      attr(result$sample[[i]], "thin") = attr(result$sample[[i]], "thin")
      attr(result$sample[[i]], "accept") = attr(result$sample[[i]], "accept")
    }
  }

  return(result)
}


#' @export
plot.bayfit = function(x, ...) plot(x$sample, ...)

#' @export
summary.bayfit = function(x, ...) summary(x$sample, ...)

#' @export
print.bayfit = function(x, ...) print(summary(x$sample, ...))

#' @export
coef.bayfit = function(x, ...) {
  sum = summary(x$sample, ...)
  sum$statistics[,1]
}

#' @export
confint.bayfit = function(x, level = 0.95, ...) {
  alpha = 1 - level
  sum = summary(x$sample, quantile = c(alpha/2, 1 - alpha/2))
  sum$quantiles
}

