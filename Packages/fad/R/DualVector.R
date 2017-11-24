
# Class definition and constructors -------------------------------------------------------------------------------

# A DualVector stores a vector of values and matrix of gradients (or Jacobian matrix)
# In order to be a proper replacement for numeric, it inherits from numeric
DualVector = setClass("DualVector", slots = list(val = "numeric", grad = "matrix"))

setMethod("initialize", 'DualVector', function(.Object, val, grad) {
  if(!inherits(grad, "matrix")) stop("grad must be a matrix")
  if(nrow(grad) != length(val)) stop("Length of val and nrow of grad should concide.")
  .Object@val = val
  .Object@grad = grad
  .Object
})

#' @export
dualvector = function(val = 0, grad = matrix(rep(0, length(val)), nrow = length(val))) new("DualVector", val, grad)

# Conversions -----------------------------------------------------------------------------------------------------

# Convert to a DualVector
#' @export
setGeneric("as.DualVector", def = function(x) {stop("Cannot convert to DualVector")})
#' @export
setMethod("as.DualVector", c("numeric"), function(x) {
  oldAttr = attributes(x)
  x = dualvector(x)
  attributes(x) = c(attributes(x), oldAttr)
  x
})
#' @export
setMethod("as.DualVector", c("logical"), function(x) {
  oldAttr = attributes(x)
  x = dualvector(as.numeric(x))
  attributes(x) = c(attributes(x), oldAttr)
  x
  })
#' @export
setMethod("as.DualVector", c("DualVector"), function(x) x)

# Concatenate DualVectors
setMethod("c", "DualVector", function(x, ..., recursive = FALSE) {
  args = list(x, ...)
  # Check that grads have all the same number of columns
  ncols = sapply(args, function(x) ncol(x@grad))
  if(length(unique(ncols)) != 1) stop("Attemp to concatenate DualVectors with different ncols for grad")
  # Check if at least one vector contains names and if so
  if(any(unlist(lapply(args, function(x) !is.null(attr(x, "names")))))) {
    names = unlist(lapply(args, function(x) {
      names = attr(x, "names")
      if(is.null(names)) names = rep("", length(x))
      names
    }))
  } else {
    names = NULL
  }
  # Stack values and matrices (by rows)
  x@val = unlist(lapply(args, function(x) x@val))
  x@grad = matrix(unlist(lapply(args, function(x) x@grad)), ncol = ncols, byrow = TRUE)
  if(!is.null(names)) attr(x, "names") = names
  x
})

# Convert to list of DualVector
#' @export
as.list.DualVector = function(x, ...) {
    out = vector("list", length(x))
    names(out) = names(x)
    for(i in 1:length(out)) out[[i]] = x[[i]]
    out
}
# setMethod("as.list", "DualVector", function(x) {
#   out = vector("list", length(x))
#   names(out) = names(x)
#   for(i in 1:length(out)) out[[i]] = x[[i]]
#   out
# })

# Test properties of a DualVector ---------------------------------------------------------------------------------

#' @export
is.DualVector = function(x) inherits(x, "DualVector")
#' @export
length.DualVector = function(x) length(x@val)
#' @export
is.finite.DualVector = function(x) is.finite(x@val)
#' @export
is.infinite.DualVector = function(x) is.infinite(x@val)
#' @export
is.na.DualVector = function(x) is.na(x@val)
#' @export
is.nan.DualVector = function(x) is.nan(x@val)

# Setters and getters ---------------------------------------------------------------------------------------------

# Retrieve value component
#' @export
setGeneric("val", function(x) {stop("Method not implemented")})
#' @export
setMethod("val", "DualVector", function(x) x@val)

# Retrieve grad component
#' @export
setGeneric("grad", function(x) {stop("Method not implemented")})
#' @export
setMethod("grad", "DualVector", function(x) x@grad)

# Extract values
#' @export
setMethod("[", "DualVector", function(x, i) {
  if(max(i) > length(x@val)) stop("Subscript out of bounds")
  x@val = x@val[i]
  x@grad = x@grad[i,,drop = FALSE]
  x
})
#' @export
setMethod("[[", "DualVector", function(x, i) {
  if(length(i) > 1) stop("Multiple indices attempted when using [[]]")
  x[i]
})

# Assign values
#' @export
setMethod("[<-", "DualVector", function(x, i, value) {
  if(is.character(i)) i = which(names(x) %in% i)
  if(max(i) > length(x@val)) stop("Subscript out of bounds")
  # Coerce to DualVector (allows assigning literal constants)
  if(!is.DualVector(value)) value = as.DualVector(value)
  x@val[i] = value@val
  x@grad[i,] = value@grad
  x
})
#' @export
setMethod("[[<-", "DualVector", function(x, i, value) {
  if(length(i) > 1) stop("Multiple indices attempted when using [[]]")
  x[i] = value
  x
})

#' @export
setMethod("names<-", c("DualVector", "ANY"), function(x, value) {
  attr(x, "names") = value
  x
})



# Print a Dualvector ----------------------------------------------------------------------------------------------

setMethod("show", "DualVector", function(object) {
  if(is.null(attr(object, "names"))) {
    for(i in 1:length(object@val))
      cat(paste0("[", i, "]"), "val: ", object@val[i], "grad:", object@grad[i,], '\n',sep = " ")
  } else {
    names = attr(object, "names")
    for(i in 1:length(object@val))
      cat(names[i], "=", "val:", object@val[i], "grad:", object@grad[i,], '\n',sep = " ")
  }
})

