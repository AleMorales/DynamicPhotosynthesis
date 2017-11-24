#' @include generics.R
NULL

# Function that calculates non-linear fit of a model. Automatically dispatched to method based on formula or black-box
#' @export
nlfit = function(model,  ...) {
  if(inherits(model, "function"))
    nlfit_bbox(model, ...)
  else
    nlfit_formula(model,...)
}

# Default print of the object
#' @export
#' @method print nlfit
print.nlfit = function(x, ...) {
  cat("Coefficients:\n")
  print(coef(x))
  cat("NLL: ")
  cat(x$NLL, "\n")
  cat("RMSE: ")
  cat(rmse(x))
}

# Extract coefficients
#' @export
#' @method coef nlfit
coef.nlfit = function(object, ...) {
  object$par
}

# Print diagnostics of the optimization
#' @export
#' @method diagnostics nlfit
diagnostics.nlfit = function(x) {
  print(x[c("algorithm", "status", "message", "nevals")])
}

# calculates residuals (observations - predictions)
#' @export
#' @method residuals nlfit
residuals.nlfit = function(object, ...) {
  object$obs - object$fitted
}

#' Calculate root mean square error of a model
#' @export
rmse = function(x) {
  sd(residuals(x), na.rm = TRUE)
}

#' Calculate mean absolute error of a model
#' @export
mae = function(x) {
  mean(abs(residuals(x)), na.rm = TRUE)
}

#' Calculate coefficient of determination of a model
#' @export
cd = function(x) {
  1 - var(residuals(x), na.rm = TRUE)/var(x$obs, na.rm = TRUE)
}

# Default print of the object
#' @export
#' @method plot nlfit
plot.nlfit = function(object, type = "1:1", var = NULL, ...) {
  if(type == "1:1") {
    plot(object$obs, object$fitted, ...)
    abline(a = 0, b = 1)
    abline(lm(object$fitted~object$obs), lty = 2)
    legend("topleft", c("1:1", "Regression"), lty = 1:2, bty = "n")
  } else if(type == "residual") {
    plot(density(residuals(object)), ...)
  } else if(type == "response") {
    plot(var, object$obs, ...)
    lines(var, object$fitted, ...)
  }
}


# Print diagnostics of the optimization
#' @export
#' @method diagnostics nlfit
diagnostics.nlfit = diagnostics.mlefit


# Calculate the likelihood profiles for different parameters
#' @export
#' @method profile nlfit
profile.nlfit = profile.mlefit

# Calculate the likelihood profiles for different parameters
#' @export
#' @method slice nlfit
slice.nlfit = slice.mlefit

# Calculate likelihood slices for different parameters
#' @export
#' @method slice nlfit
slice.nlfit = slice.mlefit

# Calculate confidence interval from a profile
#' @export
#' @method confint nlfit
confint.nlfit = confint.mlefit


# Calculate the variance covariance matrix of the model
#' @export
#' @method vcov nlfit
vcov.nlfit = vcov.mlefit
