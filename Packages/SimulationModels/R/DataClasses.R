# Inputs for simulations are stored in R6 classes. They ensure that inputs are in the correct order,
# provide with default values, perform unit conversions, etc

###################################################################################################
###################################################################################################
#################################### Setters and Getters ##########################################
###################################################################################################
###################################################################################################
# Auxilliary function to test if a name specified by the user is present in the list of inputs
check_names = function(names, vector, element) {
  present = names %in% names(vector)
  if(!all(present)) stop("There is no ", element, " called ", names[which(!present)])
}

# Auxilliary function to test if a name specified by the user is present in the list of inputs
# when the names are in the columns of the input
check_colnames = function(names, vector, element) {
  present = names %in% colnames(vector)
  if(!all(present)) stop(qq("There is no @{element} called \"@{names[which(!present)]}\"\n"))
}

# Setter to change the value of a forcing
set_forcings = function(names, values) {
  check_names(names, self$Forcings$Values,"forcing")
  self$Forcings$Values[[names]] = cbind(values[,1], values[,2]*self$Forcings$Coefs[names])
}

# Getter to get current values for a forcing
get_forcings = function(names) {
  check_names(names, self$Forcings$Values,"forcing")
  cbind(self$Forcings$Values[[names]][,1], self$Forcings$Values[[names]][,2]/self$Forcings$Coefs[names])
}

# Setter to change the value of a parameter
set_parameters = function(names, values) {
  check_names(names, self$Parameters$Values,"parameter")
  self$Parameters$Values[names] = values*self$Parameters$Coefs[names]
}

# Getter to obtain the value of a parameter
get_parameters = function(names) {
  check_names(names, self$Parameters$Values,"parameter")
  self$Parameters$Values[names]/self$Parameters$Coefs[names]
}

# Setter to change the value of a colum in the matrix of parameters associated to a parameter scan object
set_parameters_scan = function(names, values) {
  check_colnames(names, self$Parameters$Values,"parameter")
  correction = matrix(self$Parameters$Coefs[names], nrow = nrow(values), ncol = ncol(values), byrow = TRUE)
  self$Parameters$Values[,names] = values*correction
}

# Getter to obtain the value of a colum in the matrix of parameters associated to a parameter scan object
get_parameters_scan = function(names) {
  check_colnames(names, self$Parameters$Values,"parameter")
  correction = matrix(self$Parameters$Coefs[names], nrow = nrow(self$Parameters$Values), ncol = length(names), byrow = TRUE)
  self$Parameters$Values[,names]/correction
}

# Setter to change the initial value of a state
set_states = function(names, values) {
  check_names(names, self$States$Values,"state")
  self$States$Values[names] = values*self$States$Coefs[names]
}

# Getter to obtain the initial value of a state
get_states = function(names) {
  check_names(names, self$States$Values,"state")
  self$States$Values[names]/self$States$Coefs[names]
}

# Setter to change the timepoints for the next simulation
set_time = function(values) {self$Time = as.numeric(values)}

# Setter to change the settings for the next simulation
set_settings = function(names, values) {
  check_names(names, self$Settings,"setting")
  if(length(names) > 1) {
    for(i in 1:length(names)) {
      self$Settings[[names[i]]] = values[i]
    }
  } else {
    self$Settings[[names]] = values
  }

}

###################################################################################################
###################################################################################################
####################################### Check Input ###############################################
###################################################################################################
###################################################################################################

# Check that forcings passed to the constructor are correct
check_forcings = function(Forcings) {
  stopifnot(
    inherits(Forcings, "list"),
    length(Forcings) == 2,
    all(names(Forcings) %in% c("Coefs","Values")),
    inherits(Forcings[["Values"]], c("list","NULL")),
    inherits(Forcings[["Coefs"]], c("numeric","NULL")),
    if(is.null(Forcings[["Values"]])) is.null(Forcings[["Coefs"]]) else TRUE,
    if(is.null(Forcings[["Coefs"]])) is.null(Forcings[["Values"]]) else TRUE,
    length(Forcings[["Values"]]) == length(Forcings[["Coefs"]]),
    all(names(Forcings[["Values"]]) %in% names(Forcings[["Coefs"]])),
    all(plyr::laply(Forcings[["Values"]], function(x) ncol(x) == 2)),
    all(plyr::laply(Forcings[["Values"]], function(x) !is.unsorted(x[,1]))))
}

# Check that states passed to the constructor are correct
check_states = function(States) {
  stopifnot(
    inherits(States, "list"),
    length(States) == 2,
    all(names(States) %in% c("Coefs","Values")),
    inherits(States[["Values"]], "numeric"),
    inherits(States[["Coefs"]], "numeric"),
    all(names(States[["Values"]]) %in% names(States[["Coefs"]])),
    length(States[["Values"]]) == length(States[["Coefs"]]))
}

# Check that parameters passed to the constructor are correct
check_parameters = function(Parameters) {
  stopifnot(
    inherits(Parameters, "list"),
    length(Parameters) == 2,
    all(names(Parameters) %in% c("Coefs","Values")),
    inherits(Parameters[["Values"]], c("numeric","NULL")),
    inherits(Parameters[["Coefs"]], c("numeric","NULL")),
    class(Parameters[["Values"]]) == class(Parameters[["Coefs"]]),
    all(names(Parameters[["Values"]]) %in% names(Parameters[["Coefs"]])),
    length(Parameters[["Values"]]) == length(Parameters[["Coefs"]]))
}


# Check that parameters to perform a parameter scan passed to the constructor are correct
# In this case we must have a matrix of parameters. It does not make sense to have NULLs
# The names of the parameters are associated to the columns
check_parameters_scan = function(Parameters) {
  stopifnot(
    inherits(Parameters, "list"),
    length(Parameters) == 2,
    all(names(Parameters) %in% c("Coefs","Values")),
    inherits(Parameters[["Values"]], "matrix"),
    inherits(Parameters[["Coefs"]], "numeric"),
    all(colnames(Parameters[["Values"]]) %in% names(Parameters[["Coefs"]])),
    ncol(Parameters[["Values"]]) == length(Parameters[["Coefs"]]))
}


# Check that observed passed to the constructor are correct
check_observed = function(Observed) {
  stopifnot(
    inherits(Observed, "list"),
    length(Observed) == 2,
    all(names(Observed) %in% c("Coefs","Names")),
    inherits(Observed[["Names"]], c("character","NULL")),
    inherits(Observed[["Coefs"]], c("numeric","NULL")),
    if(is.null(Observed[["Names"]])) is.null(Observed[["Coefs"]]) else TRUE,
    if(is.null(Observed[["Coefs"]])) is.null(Observed[["Names"]]) else TRUE,
    all(Observed[["Names"]] %in% names(Observed[["Coefs"]])),
    length(Observed[["Names"]]) == length(Observed[["Coefs"]]))
}

# Check that the time passed to constructor are correct
check_time = function(Time) {
  stopifnot(
    is.numeric(Time),
    !is.unsorted(Time),
    length(Time) > 1
  )
}

# Creates Settings by populating it with default values when not specified
# This allows creating correct Settings whenever a list or NULL is passed
create_settings_ODE = function(Settings, nstates, nobs, nparm) {
  stopifnot(inherits(Settings, c("list","NULL")))
  if(is.null(Settings[["atol"]])) Settings[["atol"]] <- 1e-6
  if(is.null(Settings[["rtol"]])) Settings[['rtol']] <- 1e-6
  if(is.null(Settings[["maxsteps"]])) Settings[["maxsteps"]] <- 500
  if(is.null(Settings[["method"]])) Settings[["method"]] <- "bdf"
  if(Settings[["method"]] == "bdf") {
    if(is.null(Settings[["maxord"]])) Settings[["maxord"]] <- 5
  } else {
    if(is.null(Settings[["maxord"]])) Settings[["maxord"]] <- 12
  }
  if(is.null(Settings[["hini"]])) Settings[["hini"]] <- 0
  if(is.null(Settings[["hmin"]])) Settings[["hmin"]] <- 0
  if(is.null(Settings[["hmax"]])) Settings[["hmax"]] <- 0
  if(is.null(Settings[["maxerr"]])) Settings[["maxerr"]] <- 7
  if(is.null(Settings[["maxnonlin"]])) Settings[["maxnonlin"]] <- 3
  if(is.null(Settings[["maxconvfail"]])) Settings[["maxconvfail"]] <- 10
  if(is.null(Settings[["maxtime"]])) Settings[["maxtime"]] <- 60
  if(is.null(Settings[["which_states"]])) Settings[["which_states"]] <- 1:nstates
  if(is.null(Settings[["which_observed"]])) Settings[["which_observed"]] <- if(nobs > 0) 1:nobs else vector("integer")
  if(is.null(Settings[["positive"]])) Settings[["positive"]] = 1
  if(is.null(Settings[["minimum"]])) Settings[["minimum"]] = -1e-6
  if(is.null(Settings[["stability"]])) Settings[["stability"]] = T
  if(is.null(Settings[["sensitivity"]])) Settings[["sensitivity"]] = F
  if(is.null(Settings[["sens_type"]])) Settings[["sens_type"]] = "simultaneous"
  if(is.null(Settings[["sens_fun"]])) Settings[["sens_fun"]] = F
  if(is.null(Settings[["sens_error"]])) Settings[["sens_error"]] = F
  if(is.null(Settings[["which_sens"]])) Settings[["which_sens"]] = 1:nparm
  if(is.null(Settings[["reset"]])) Settings[["reset"]] = vector("integer")
  if(is.null(Settings[["silent"]])) Settings[["silent"]] = FALSE
  if(is.null(Settings[["number_resets"]])) Settings[["number_resets"]] = 2L
  if(is.null(Settings[["force_positive"]])) Settings[["force_positive"]] = FALSE
  if(is.null(Settings[["which_positive"]])) Settings[["which_positive"]] = as.integer(1:nstates)
  Settings[["jacobian"]] = 0
  return(Settings)
}

create_settings_DAE = function(Settings, nstates, nobs) {
  stopifnot(inherits(Settings, c("list","NULL")))
  #if(is.null(Settings[["vars_id"]])) stop("You must provide vars_id in settings")
  if(is.null(Settings[["atol"]])) Settings[["atol"]] <- 1e-6
  if(is.null(Settings[["rtol"]])) Settings[['rtol']] <- 1e-6
  #if(is.null(Settings[["suppress_alg"]])) Settings[["suppress_alg"]] = 1
  if(is.null(Settings[["calcic"]])) Settings[["calcic"]] = 1
  if(is.null(Settings[["maxsteps"]])) Settings[["maxsteps"]] <- 500
  if(is.null(Settings[["maxord"]])) Settings[["maxord"]] <- 5
  if(is.null(Settings[["hini"]])) Settings[["hini"]] <- 0
  if(is.null(Settings[["hmax"]])) Settings[["hmax"]] <- 0
  if(is.null(Settings[["maxerr"]])) Settings[["maxerr"]] <- 7
  if(is.null(Settings[["maxnonlin"]])) Settings[["maxnonlin"]] <- 3
  if(is.null(Settings[["maxconvfail"]])) Settings[["maxconvfail"]] <- 10
  if(is.null(Settings[["which_states"]])) Settings[["which_states"]] <- 1:nstates
  if(is.null(Settings[["which_observed"]])) Settings[["which_observed"]] <- if(nobs > 0) 1:nobs else vector("integer")
  if(is.null(Settings[["positive"]])) Settings[["positive"]] = 1
  if(is.null(Settings[["minimum"]])) Settings[["minimum"]] = -1e-6
  Settings[["jacobian"]] = 0
  return(Settings)
}

# Check that we have a valid, non-null pointer
check_pointer = function(pointer) {
  null_pointer = new("externalptr")
  attributes(null_pointer) = attributes(pointer)
  stopifnot(
    inherits(pointer, "NativeSymbol"),
    storage.mode(pointer) == "externalptr",
    !identical(pointer, null_pointer))
}

###################################################################################################
###################################################################################################
####################################### Constructors ##############################################
###################################################################################################
###################################################################################################

# We use the constructor to check validity.
initialize_ODE = function(Forcings = list(Values = NULL, Coefs = NULL),
                         States = list(Values = NULL, Coefs = NULL),
                         Parameters = list(Values = NULL, Coefs = NULL),
                         Observed = list(Names = NULL, Coefs = NULL),
                         Time = 1:2,
                         Settings = NULL,
                         model_function = NULL) {
  # Check that forcings are correct and assign them
  check_forcings(Forcings); self$Forcings = Forcings
  check_states(States); self$States = States
  check_parameters(Parameters); self$Parameters = Parameters
  check_observed(Observed); self$Observed = Observed
  check_time(Time); self$Time = as.numeric(Time)
  Settings = create_settings_ODE(Settings, length(States$Values), length(Observed$Names), length(Parameters$Values)); self$Settings = Settings
  check_pointer(eval(parse(text = model_function))); self$model_function = parse(text = model_function)
}


initialize_DAE = function(Forcings = list(Values = NULL, Coefs = NULL),
                          States = list(Values = NULL, Coefs = NULL),
                          Derivatives = list(Values = NULL, Coefs = NULL),
                          Parameters = list(Values = NULL, Coefs = NULL),
                          Observed = list(names = NULL, Coefs = NULL),
                          Time = 1:2,
                          Settings = NULL,
                          model_function = NULL) {
  check_forcings(Forcings); self$Forcings = Forcings
  check_states(States); self$States = States
  check_states(Derivatives); self$Derivatives = Derivatives
  check_parameters(Parameters); self$Parameters = Parameters
  check_observed(Observed); self$Observed = Observed
  check_time(Time); self$Time = as.numeric(Time)
  Settings = create_settings_DAE(Settings, length(States$Values), length(Observed$Names)); self$Settings = Settings
  check_pointer(eval(parse(text = model_function))); self$model_function = parse(text = model_function)
}

###################################################################################################
###################################################################################################
##################################### copy function ###########################################
###################################################################################################
###################################################################################################
copy <- function(x, ...) UseMethod("copy", x)

copy.DAEmodel = function(old) old$clone() #{
#   DAEmodel$new(
#    States = old$States,
#    Derivatives = old$Derivatives,
#    Parameters = old$Parameters,
#    Forcings = old$Forcings,
#    Observed = old$Observed,
#    Time = old$Time,
#    Settings = old$Settings,
#    model_function = old$model_function,
#   )
# }

 copy.ODEmodel = function(old) old$clone() #{

###################################################################################################
###################################################################################################
##################################### Class definitions ###########################################
###################################################################################################
###################################################################################################


# Class to store data and methods associated to ODE models
#' @export DAEmodel
DAEmodel = R6::R6Class("DAEmodel",
                   public = list(
                     # Inputs of the model
                     Forcings = NA,
                     States = NA,
                     Derivatives = NA,
                     Parameters = NA,
                     Observed = NA,
                     Time = NA,
                     Settings = NA,
                     # Function that define a model
                     model_function = NA,
                     # Method to initialize an instance
                     initialize = initialize_DAE,
                     # Assign values to the inputs of the model with the proper unit conversion coefficients
                     set_forcings = set_forcings,
                     get_forcings = get_forcings,
                     set_parameters = set_parameters,
                     get_parameters = get_parameters,
                     set_states = set_states,
                     get_states = get_states,
                     set_time = set_time,
                     set_settings = set_settings
                   ))


# Class to store data and methods associated to DAE models
#' @export ODEmodel
ODEmodel = R6::R6Class("ODEmodel",
                   inherit = DAEmodel,
                   public = list(
                     # Method to initialize an instance
                     initialize = initialize_ODE
                   ))
