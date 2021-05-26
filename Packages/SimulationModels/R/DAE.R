
###################################################################################################
###################################################################################################
##################################### DAEmodel methods ############################################
###################################################################################################
###################################################################################################

# Run a dynamic simulation
#' @export
ida.DAEmodel = function(model) {
    sim = ida_Cpp_stl(times = model$Time, states = model$States$Values,
                derivatives = model$Derivatives$Values,
                parameters = model$Parameters$Values,
                forcings_data = model$Forcings$Values,
                settings = model$Settings,
                model = eval(model$model_function),
                jacobian = eval(model$model_function))
  colnames(sim) = c("time", names(model$States$Values)[Settings[["which_states"]]],
                    model$Observed$Names[Settings[["which_observed"]]])
  corrections = matrix(c(model$States$Coefs[Settings[["which_states"]]],
                         model$Observed$Coefs[Settings[["which_observed"]]]),
                       nrow = nrow(sim), ncol = length(Settings[["which_states"]]) +
                         length(Settings[["which_observed"]]), byrow = TRUE)
  sim[,-1] = sim[,-1]/corrections
  class(sim) = c("OutputMatrix")
  sim
}

# Need to define steady

# Need to define analyse_steady
