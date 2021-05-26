# This will implements tests in the future
test = function() {print("Tests are not implemented yet")}

###################################################################################################
###################################################################################################
########################################## Plots ##################################################
###################################################################################################
###################################################################################################

# Matrix of outputs
OutputMatrix <- structure(matrix, class = "OutputMatrix")

# Each variable vs time in its own plot
# Make sure the names are maintained there
#' @export
plot.OutputMatrix = function(mat, subset = 1:nrow(mat), which = 2:ncol(mat), mfrow = c(3,3)) {
  par(mfrow = mfrow, las = 1)
  for(i in which) {
    ylab = if(inherits(i,"character")) {i} else {colnames(mat)[i]}
    plot(mat[subset,1], mat[subset,i], xlab = "Time", ylab = ylab, type = "l")
  }
}

# Array of outputs
OutputArray <- structure(array, class = "OutputArray")

# Each variable vs time in its own plot but multiple replicates with colors
# Make sure the names are maintained there
# Add legend for first one
#' @export
plot.OutputArray = function(arr, subset = 1:nrow(arr), which = 2:ncol(arr), mfrow = c(3,3)) {
  if(length(dim(arr)) > 3) stop("I cannot plot arrays with 4 or more dimensions.")
  par(mfrow = mfrow, las = 1)
  for(i in which) {
    ylab = if(inherits(i,"character")) {i} else {colnames(arr)[i]}
    ylim = c(0.9,1.1)*range(arr[,i,])
    for(j in 1:dim(arr)[3]) {
      if(j == 1) {
        plot(arr[subset,1,j], arr[subset,i,j], xlab = "Time", ylab = ylab, type = "l", ylim = ylim)
      } else {
        lines(arr[subset,1,j], arr[subset,i,j], col = j)
      }
    }
  }
}

# Generic functions to run simulations
#' @export
cvode <- function(x, ...) UseMethod("cvode", x)

#' @export
ida <- function(x, ...) UseMethod("ida", x)


###################################################################################################
###################################################################################################
########################################## Convert ################################################
###################################################################################################
###################################################################################################

#' @export
as.data.frame.OutputMatrix = function(x) {
  class(x) = "matrix"
  as.data.frame(x)
}

#' @export
as.matrix.OutputMatrix = function(x) {
  class(x) = "matrix"
  x
}






























