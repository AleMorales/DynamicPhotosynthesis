plot.Simulation = function(output, ...) {
  class(output) = "deSolve"
  deSolve:::plot.deSolve(output)
}


