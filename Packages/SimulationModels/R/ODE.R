###################################################################################################
###################################################################################################
##################################### ODEmodel methods ############################################
###################################################################################################
###################################################################################################

# Run a dynamic simulation
cvode.ODEmodel = function(model) {
    sim = wrap_cvodes(times = model$Time, states = model$States$Values,
                parameters = model$Parameters$Values,
                forcings_data = model$Forcings$Values,
                settings = model$Settings,
                model = eval(model$model_function),
                jacobian = eval(model$model_function))

  if(model$Settings[["sensitivity"]]) {
    NS = length(model$Settings[["which_sens"]])
    NEQ = length(model$Settings[["which_states"]])
    names_SC = paste0("S_",rep(names(model$States$Values)[model$Settings[["which_states"]]], each = NS),"_",
                      rep(names(model$Parameters$Values)[model$Settings[["which_sens"]]], NEQ))

    colnames(sim) = c("time", names(model$States$Values)[model$Settings[["which_states"]]],
                      model$Observed$Names[model$Settings[["which_observed"]]], names_SC)
    corrections = matrix(c(model$States$Coefs[model$Settings[["which_states"]]],
                           model$Observed$Coefs[model$Settings[["which_observed"]]],
                           rep(model$States$Coefs[model$Settings[["which_states"]]], each = NS)/
                           rep(model$Parameters$Coefs[model$Settings[["which_sens"]]], NEQ)),
                         nrow = nrow(sim), ncol = length(model$Settings[["which_states"]])*(1 + NS) +
                           length(model$Settings[["which_observed"]]), byrow = TRUE)

  } else {
    colnames(sim) = c("time", names(model$States$Values)[model$Settings[["which_states"]]],
                      model$Observed$Names[model$Settings[["which_observed"]]])
    corrections = matrix(c(model$States$Coefs[model$Settings[["which_states"]]],
                           model$Observed$Coefs[model$Settings[["which_observed"]]]),
                         nrow = nrow(sim), ncol = length(model$Settings[["which_states"]]) +
                           length(model$Settings[["which_observed"]]), byrow = TRUE)
  }


  sim[,-1] = sim[,-1]/corrections
  class(sim) = c("OutputMatrix")
  sim
}

# Need to define steady

# Need to define analyse_steady
