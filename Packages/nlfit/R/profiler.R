#' @importFrom foreach foreach %do%
NULL

profiler = function(object, parm = names(object$par), options = object$options, algorithm = object$algorithm,
                    lower = object$lower, upper = object$upper,
                    level = 0.95, lower_delta = 0.2, upper_delta = lower_delta, nstep = 10, trace = FALSE,
                    direction = "both", maxsteps = 10, tolerance = 0, cl = NULL,
                    search_bounds = TRUE, eval = FALSE, ...) {

  # For one-parameter models use a slice
  if(length(object$par) == 1) eval = TRUE


  # Make sure that all names in parm are names in object$par
  if(!all(parm %in% names(object$par)))
    stop("Profiles can only be calculated for fitted model parameters as reported by the coef() method")

  # Allows to specify partial changes to options
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

  # Reference values of NLL and fitted parameters
  NLL = object$NLL
  npar = object$par[parm]

  # Indices for parameters for which a profile is being calculated
  ipar = numeric(length(parm))
  for(i in 1:length(parm))
    ipar[i] = which(names(object$par) == parm[i])

  # Difference in NLL required for desired confidence level
  dNLL = qchisq(p = level, df = 1)/2
  maxdNLL = 5*dNLL # To detect excessive initial delta

  # Define a delta for each parameter and each half of the profile
  # the delta vectors could be named or not. If not named, it must be 1 or length(parm)
  # if named, it must only contain names within parm
  if(is.null(names(lower_delta))) {
    if(length(lower_delta) == 1)
      lower_delta = rep(lower_delta, times = length(parm))
    else if(length(lower_delta) != length(parm))
      stop("Length of lower_delta and parm do not coincide.")
  } else {
    if(!all(names(lower_delta) %in% parm))
      stop("Argument lower_delta can only specify values for parameters within argument parm")
    input_lower_delta = lower_delta
    lower_delta = rep(0.2, times = length(parm))
    names(lower_delta) = parm
    lower_delta[names(input_lower_delta)] = unname(input_lower_delta)
  }

  if(is.null(names(upper_delta))) {
    if(length(upper_delta) == 1)
      upper_delta = rep(upper_delta, times = length(parm))
    else if(length(upper_delta) != length(parm))
      stop("Length of upper_delta and parm do not coincide.")
  } else {
    if(!all(names(upper_delta) %in% parm))
      stop("Argument upper_delta can only specify values for parameters within argument parm")
    input_upper_delta = upper_delta
    upper_delta = rep(0.2, times = length(parm))
    names(upper_delta) = parm
    upper_delta[names(input_upper_delta)] = unname(input_upper_delta)
  }

  # Generate a likelihood profile that includes dNLL on each side, for parm
  # If a cluster is provided use it to calculate the profiles in parallel
  calc_profile = function(i, ...) {
    if(exists("lNLL")) rm(lNLL)
    if(exists("iNLL")) rm(iNLL)
    if(exists("uNLL")) rm(uNLL)

    # Adjust start value and fixed the parameter used for profiling
    start  = object$all_par

    # Taken into account that lower and upper may not contain this particular parameter
    low = max(valid_lower[parm[i]], -Inf, na.rm = T)
    up = min(valid_upper[parm[i]], Inf, na.rm = T)

    # Ensure that the lowest position contains a difference of at least dNLL but does not exceed maxdNLL
    # Use a bisection method to fall within the range
    # The exception  is when initial value is below dNLL. In that case double the delta
    if(direction  %in% c("both", "lower")) {
      lpos = max((1 - lower_delta[i])*object$par[parm[i]], low)
      n = 0
      last_pos = object$par[parm[i]]
      while(TRUE) {
        if(n > maxsteps) {
          warning("Maximum number of steps reached when calculating boundary for lower half-profile")
          break
        }
        n = n + 1
        if(exists("lNLL")) start = lNLL$all_par else start= object$all_par
        start[parm[i]] = lpos
        lNLL = update(object, start = start, options = object$options, algorithm = algorithm,
                      lower = valid_lower, upper = valid_upper,
                      fixed  = c(object$fixed, parm[i]), trace = trace, eval = eval, ...)
        if(is.na(lNLL$NLL)) {
          stop("Error when calculating lower boundary")
        }
        if((lNLL$NLL < NLL) & (abs(lNLL$NLL - NLL) > tolerance)) {
          stop("A better fit was obtained for par = ", paste(names(lNLL$all_par), " = ", lNLL$all_par, collapse = " "), ". Please rerun optimization.")
        }
        if(search_bounds & (lNLL$NLL - NLL) > maxdNLL & lpos > low) {
          if(last_pos > lpos)
            new_lpos = (lpos + last_pos)/2
          else
            new_lpos = (lpos + object$par[parm[i]])/2
          last_pos = lpos
          lpos = new_lpos
        } else if(search_bounds & (lNLL$NLL - NLL) < dNLL & lpos > low) {
          if(last_pos < lpos)
            new_lpos = (lpos + last_pos)/2
          else
            new_lpos = object$par[parm[i]] - 2*(object$par[parm[i]] - lpos)
          last_pos = lpos
          lpos = new_lpos
        } else if ((lNLL$NLL - NLL) < dNLL & lpos <= low)
          stop("Error when calculating profile for parameter ", parm[i], ". The lower boundary is within the confidence interval requested.")
        else
          break
      }

    }

    if(direction %in% c("both", "upper")) {
      # Ensure that the highest position contains a difference of at least dNLL
      upos = min((1 + upper_delta[i])*object$par[parm[i]], up)
      n = 0
      last_pos = object$par[parm[i]]
      while(TRUE) {
        if(n > maxsteps) {
          warning("Maximum number of steps reached when calculating boundary for upper half-profile")
          break
        }
        n = n + 1
        if(exists("uNLL")) start = uNLL$all_par else start= object$all_par
        start[parm[i]] = upos
        uNLL = update(object, start = start, options = object$options, algorithm = algorithm,
                      lower = valid_lower, upper = valid_upper,
                      fixed  = c(object$fixed, parm[i]), trace = trace, eva = eval, ...)
        if(is.na(uNLL$NLL)) {
          stop("Error when calculating upper boundary")
        }
        if((uNLL$NLL < NLL) & (abs(uNLL$NLL - NLL) > tolerance))
          stop("A better fit was obtained for par = ", paste(names(uNLL$all_par), " = ", uNLL$all_par, collapse = " "), ". Please rerun optimization.")
        if(search_bounds & (uNLL$NLL - NLL) > maxdNLL & upos < up) {
          if(last_pos < upos)
            new_upos = (upos + last_pos)/2
          else
            new_upos = (upos + object$par[parm[i]])/2
          last_pos = upos
          upos = new_upos
        } else if(search_bounds & (uNLL$NLL - NLL) < dNLL & upos < up) {
          if(last_pos > upos)
            new_upos = (upos + last_pos)/2
          else
            new_upos = object$par[parm[i]] + 2*(upos - object$par[parm[i]])
          last_pos = upos
          upos = new_upos
        } else if ((uNLL$NLL - NLL) < dNLL & upos >= up)
          stop("Error when calculating profile for parameter ", parm[i], ". The upper boundary is within the confidence interval requested.")
        else
          break
      }
    }

    # Calculate 2 half-profiles of length nsteps with upper and lower positions determined by lpos and upos
    profile = matrix(NA, nrow = nstep, ncol = 2)
    colnames(profile) = c("par", "DNLL")
    if(direction %in% c("both", "lower")) {
      init = 1
      profile[1:(nstep/2),1] = seq(object$par[parm[i]], lpos,length.out = nstep/2 + 1)[-1]
      profile[1,2] = lNLL$NLL - NLL
    } else {
      init = nstep/2 + 1
    }
    if(direction %in% c("both", "upper")) {
      end = nstep
      profile[(nstep/2 + 1):nstep,1] = seq(object$par[parm[i]], upos, length.out = nstep/2 + 1)[-1]
      profile[nrow(profile), 2] = uNLL$NLL - NLL
    } else {
      end = nstep/2
    }
    for(j in init:end){
      if(exists("iNLL") & (j != (nstep/2 + 1))) {
        start = iNLL$all_par
      } else {
        start = object$all_par
      }
      start[parm[i]] = profile[j,1]
      iNLL = try(update(object, start = start, options = object$options, algorithm = algorithm,
                        lower = valid_lower, upper = valid_upper,
                        fixed  = c(object$fixed, parm[i]), trace = trace, eval = eval, ...))
      if(inherits(iNLL, "try-error")) {
        profile[j,1] = profile[j,1] + 0.1*diff(profile[c(j, j + 1),1])
        start[parm[i]] = profile[j,1]
        iNLL = try(update(object, start = start, options = object$options, algorithm = algorithm,
                          lower = valid_lower, upper = valid_upper,
                          fixed  = c(object$fixed, parm[i]), trace = trace, eva = eval, ...))
      }
      profile[j,2] = iNLL$NLL - NLL
    }

    # Include nominal parameters and store profile in the list
    profile = rbind(profile[(nstep/2):1,], c(npar[parm[i]], 0), profile[(nstep/2 + 1):nstep,])

    # Calculate tail probabilities from chi-sq test (as required for calculating confidence intervals later)
    profile = cbind(profile, p = pchisq(q = profile[,2]*2, df = 1))
  }
  if(is.null(cl))
    profiles = foreach(i = 1:length(parm)) %do% calc_profile(i, ...)
  else {
    profiles = foreach(i = 1:length(parm)) %dopar% calc_profile(i, ...)
  }

  # S3 class to overload confint and plot
  names(profiles) = parm
  class(profiles) = "profile"
  attr(profiles, "level") = 0.95
  profiles

}



# Calculate the likelihood profiles for different parameters
#' @export
#' @method plot profile
plot.profile = function(x) {
  nlen = length(x)
  if(nlen == 1)
    par(mfrow = c(1,1))
  else if(nlen == 2)
    par(mfrow = c(1,2))
  else if(nlen < 5)
    par(mfrow = c(2,2))
  else
    par(mfrow = c(3,2))
  par(mgp = c(2.5,1,0))
  for(i in 1:nlen) {
      plot(x[[i]][,-3], type = "o", xlim = range(x[[i]][,1], na.rm = TRUE),
           ylim = range(x[[i]][,2], na.rm = TRUE),
           ylab = expression(Delta*log(L(theta*"|"*data))), xlab = names(x)[i])
      abline(h = qchisq(p = attr(x, "level"), df = 1)/2, lty = 2, col = 2)
  }
}

# Calculate confidence interval from a profile
#' @export
#' @method confint profile
confint.profile = function(profile, level = attr(profile, "level")) {

  n = length(profile)

  # Difference in NLL required for desired confidence level
  dNLL = qchisq(p = level, df = 1)/2

    # Setup matrix with output
  limits = as.data.frame(matrix(NA, ncol = 3, nrow = n))
  colnames(limits) = c("value", "lower", "upper")
  #limits[,"parm"] = names(profile)
  row.names(limits) = names(profile)
  limits[,"value"] = sapply(profile, function(x) x[(nrow(x) - 1)/2 + 1, 1])
  for(i in 1:n) {
    nstep = nrow(profile[[i]]) - 1
  # Find value at dNLL by interpolating from fitted cubic spline
    if(all(is.na(profile[[i]][1:(nstep/2),2]))) {
      limits[i,2] = NA
    } else {
      limits[i,2] = spline(profile[[i]][1:(nstep/2),2], profile[[i]][1:(nstep/2),1], xout = dNLL)$y
    }
    if(all(is.na(profile[[i]][(nstep/2 + 1):(nstep + 1),2]))) {
      limits[i,3] = NA
    } else {
      limits[i,3] = spline(profile[[i]][(nstep/2 + 1):(nstep + 1),2], profile[[i]][(nstep/2 + 1):(nstep + 1),1], xout = dNLL)$y
    }
  }
  limits
}
