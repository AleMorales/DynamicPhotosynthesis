# Default options, with extra settings for specific algorithms
# nlopt_algorithm = Default nlopt algorithm (may change with update_nlopt_algorithm)
# fnscale = Scalar that multiplies the output of the objective function (use -1 to switch to maximization)
# parscale = Vector that multiplies the parameters before passed to the objective function (use it to scale parameters)
# eqatol = Absolute tolerance to check equality constraints
# ineqatol = Absolute tolerance to check inequality constraints
# xatol = Absolute tolerance convergence in parameter value. Either scalar or vector of length = length(x0)
# xrtol = Relative tolerance convergence in parameter value
# fatol = Absolute tolerance convergence objective function
# frtol = Relative tolerance convergence objective function
# stopval: Minimum value of the objective function at which minimization should stop
# maxiter = Maximum number of iterations of the algorithm
# maxtime = Maximum time allowed for the algorithm to find a solution in seconds (<= 0 means no time limitation)
# dx =   Absolute initial step size for derivative-free algorithms (optional)
#        Either scalar or vector of length = length(x0)
#        An empty vector (numeric(0)) indicates no user-define initial step size
# Some parameters are not assigned by default but may be useful:
# randseed: Seed for algorithms that make use of randon numbers (if empty seed will be base on time)
# pop: Size of the population used in evolutionary algorithms (if pop = 0, internal heuristics are used)
# M: number of gradients to be remembered by truncated newton algorithm (if M = 0, internal heuristics are used)
default_options = function(algorithm) {

  options = list()

  # Common options to all algorithms
  alg = algorithms[[algorithm]]
  if(is.null(alg)) stop(paste("Algorithm", algorithm, "does not exist"))
  options$nlopt_algorithm = alg
  options$algorithm = algorithm
  options$eqatol = 1e-8
  options$ineqatol = 1e-8
  options$xatol = 0
  options$xrtol = sqrt(.Machine$double.eps)
  options$fatol = 0
  options$frtol = 0
  options$stopval = -Inf
  options$maxiter = 1e3
  options$maxtime = 0
  options$fnscale = 1
  options$parscale = 1

  # Only for derivative-free algorithms but required to be there
  if(algorithm %in% df_algorithms) options$dx = numeric(0)


  # Additional options to choose variants
  switch(algorithm,
         direct = {options$scaling = TRUE
         options$local = TRUE
         options$randomize = FALSE
         options = update_nlopt_algorithm(options)},

         stogo = {options$randomize = FALSE
         options = update_nlopt_algorithm(options)},

         var = {options$rank = 1
         options = update_nlopt_algorithm(options)},

         tnewton = {options$precondition = FALSE
         options$restart = FALSE
         options = update_nlopt_algorithm(options)}
  )

  # Additional options for specific languages
  if(algorithm == "tnewton") options$M = 0
  if(algorithm %in% c("stogo", "praxis")) options$randseed = integer()

  if(algorithm %in% c("crs", "isres", "esch")) {
    options$pop = 0
    options$randseed = integer()
  }

  return(options)
}

# If an algorithm has variants, update nlopt_algorithm based on options
update_nlopt_algorithm = function(options) {

  # Variations of the direct algorithm
  if(options$algorithm == "direct") {
    nlopt_algorithm = "GN_DIRECT"
    if(options$local) nlopt_algorithm = paste0(nlopt_algorithm, "_L")
    if(options$randomize) nlopt_algorithm = paste0(nlopt_algorithm, "_RAND")
    if(!options$scaling) nlopt_algorithm = paste0(nlopt_algorithm, "_NOSCAL")
    accepted = paste0("GN_DIRECT", c("_L_RAND", "_L", "", "_L_RAND_NOSCAL", "_L_NOSCAL", "_NOSCAL"))
    if(!(nlopt_algorithm %in% accepted))
      stop("The combination of options 'local', 'randomize' and 'scaling' for algorithm 'direct' are not correct.
           Check ?nlopt for details.")
    options$nlopt_algorithm = nlopt_algorithm
    return(options)
  }

  # Variations of the stogo algorithm
  if(options$algorithm == "stogo") {
    if(options$randomize)
      options$nlopt_algorithm = "GD_STOGO"
    else
      options$nlopt_algorithm = "GD_STOGO_RAND"
    return(options)
  }

  # Variatons of var algorithm
  if(options$algorithm == "var") {
    if(options$rank == 1) {
      options$nlopt_algorithm = "LD_VAR1"
    } else if(options$rank == 2){
      options$nlopt_algorithm = "LD_VAR2"
    } else {
      stop("The option 'rank' must be '1' or '2'. Check ?nlopt for details.")
    }
    return(options)
  }

  # Variations of tnewton algorithm
  if(options$algorithm == "tnewton") {
    nlopt_algorithm = "LD_TNEWTON"
    if(options$precondition) nlopt_algorithm = paste0(nlopt_algorithm, "_PRECOND")
    if(options$restart) nlopt_algorithm = paste0(nlopt_algorithm, "_RESTART")
    options$nlopt_algorithm = nlopt_algorithm
    return(options)
  }

  # In case algorithm does not have variants
  return(options)
}


algorithms = list(direct = "GN_DIRECT", # Has variants
                  stogo = "GD_STOGO", # Has variants
                  sqp = "LD_SLSQP",
                  lbfgs = "LD_LBFGS",
                  praxis = "LN_PRAXIS",
                  var = "LD_VAR1", # Has variants
                  tnewton = "LD_TNEWTON", # Has variants
                  crs = "GN_CRS2_LM",
                  mma = "LD_MMA",
                  cobyla = "LN_COBYLA",
                  newuoa = "LN_NEWUOA_BOUND",
                  nelder = "LN_NELDERMEAD",
                  subplex = "LN_SBPLX",
                  bobyqa = "LN_BOBYQA",
                  isres = "GN_ISRES",
                  esch = "GN_ESCH")

df_algorithms = c("direct", "praxis", "crs", "cobyla", "newuoa", "nelder", "subplex", "bobyqa", "isres", "esch")
