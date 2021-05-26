# This automatically loads all the dependencies (I used Depends)
library(SimulationModels)

###################################################################################################
###################################################################################################
########################################  Pointers  ###############################################
###################################################################################################
###################################################################################################

# Store a string to the pointer, not the pointer itself
example_model = 'getNativeSymbolInfo(name = "example_model_stl",PACKAGE = "RcppSundials")$address'

example_jacobian = 'getNativeSymbolInfo(name = "example_jacobian_stl",PACKAGE = "RcppSundials")$address'

example_dae = 'getNativeSymbolInfo(name = "example_dae_stl",PACKAGE = "RcppSundials")$address'

###################################################################################################
###################################################################################################
######################################  ODE example  ##############################################
###################################################################################################
###################################################################################################

# Create an ODE model function generator
generate_ode_example = function() {
  ODEmodel$new(
            States = list(Values = c(a = 1, b = 1, c = 1, d = 1, e = 1),
                          Coefs = c(a = 1, b = 1, c = 1, d = 1, e = 1)),
            Parameters = list(Values = c("rate"= 0.1),
                              Coefs = c("rate"= 1)),
            Forcings = list(Values = list(Forc = cbind(1:3600,1:3600)),
                            Coefs = c(Forc = 1)),
            Time = 1:5,
            Observed = list(Names = c("Forc"),
                            Coefs = c("Forc" = 1)),
            Settings = list(rtol = 1e-6,
                            atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 0,
                            hmin = 0, hmax = 0, maxerr = 5, maxnonlin = 10,
                            maxconvfail = 10, method = "bdf",jacobian = 0,
                            positive = 0, sensitivity = T),
            model = example_model)
}

# Create an instance of the model
example_ODE_model = generate_ode_example()

# Run simulation
simulation = cvode(example_ODE_model)





###################################################################################################
###################################################################################################
############################  Need to finish the rest...  #########################################
###################################################################################################
###################################################################################################







# Run analysis
example_ODE_model$set_time(seq(1,3600,by = 100))
analysis = analyse_dynamic(example_ODE_model, states = TRUE, derivatives = TRUE, jacobian = TRUE)

# Run several simulations for different set of parameters
# with an ODEscan object
example_ODE_scan = create_ODEscan(example_ODE_model, cbind(rate = log10(seq(1,2.5,length.out = 100))))
example_ODE_scan$set_time(seq(1,500,by = 5))
example_ODE_scan$set_settings("cores", 8)
multi_simulations = simulate(example_ODE_scan)
plot(multi_simulations, mfrow = c(2,3))

###################################################################################################
###################################################################################################
######################################  DAE example  ##############################################
###################################################################################################
###################################################################################################

# Create an ODE model function generator
generate_dae_example = function() {
  MiniModelDAE <- DAEmodel$new(
            States = list(Values = c(a = 1, b = 0, c = 0),
                          Coefs = c(a = 1, b = 1, c = 1)),
            Derivatives = list(Values = c(a = -0.04, b = 0.04, c = 0),
                          Coefs = c(a = 1, b = 1, c = 1)),
            Parameters = list(Values = c("rate"= 0.1),
                              Coefs = c("rate"= 1)),
            Forcings = list(Values = list(Forc = cbind(1:3600,sin((1:3600)/100*pi))),
                            Coefs = c(Forc = 1)),
            Time = 1:3600,
            Observed = list(Names = c("Forc"),
                            Coefs = c("Forc" = 1)),
            Settings = list(rtol = 1e-10,
                            atol = 1e-10, maxsteps = 1e3, maxord = 5, hini = 0,
                            hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                            maxconvfail = 10, maxtime = 0,
                            observer = 1, nder = 1, vars_id = c(1.0,1.0,0.0),
                            suppress_alg = 1, calcic = 1),
            model = example_dae)
}

# Create an instance of the model
example_DAE_model = generate_dae_example()

# Run simulation
simulation = ida(example_DAE_model)


