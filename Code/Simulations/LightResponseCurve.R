
# Load libraries ----------------------------------------------------------
library(MiniModel)
library(ggplot2)
library(readr)
library(dplyr)
library(doParallel)
library(Hmisc)

# Load data ---------------------------------------------------------------

# Parameters fitted to literature
load("Intermediate/ParametersLiterature.RData")
load("Intermediate/ParametersExperiment.RData")
parameters = params_literature
for(i in row.names(params_experiment)) {
  parameters[i,] = params_experiment[i,]
}
namesMutants = c("TPU25_spsa", "RCA_rca2", "Jmax25_rca2", "Jmax25_npq4", "Jmax25_rwt43")
parameters_mutants = parameters[namesMutants,"value"]
names(parameters_mutants) = namesMutants
parameters = parameters[-which(row.names(parameters) %in% namesMutants),]
parnames = row.names(parameters)
parameters = parameters[,"value"]
names(parameters) = parnames
parameters["gamma3"]  = 1 - parameters["gamma2"] - parameters["gamma1"]
rm(data_summary, params_experiment, params_literature, experiment_summary)


# Light Response Curve ----------------------------------------------------

# Function to simulate an induction curve given a model object
run_lrc = function(model) {
  
  # 1. Run adaptation during acclimation
  
  # Assign constant conditions during acclimation
  model$set_forcings("Ib", cbind(c(0,1), c(0,0)))
  model$set_forcings("Ir", cbind(c(0,1), c(0,0)))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("H2OR", cbind(c(0,1), c(20, 20)))
  model$set_forcings("CO2R", cbind(c(0,1), c(400, 400)))
  model$set_forcings("Ta", cbind(c(0,1), c(298.15, 298.15)))
  model$set_forcings("Tl", cbind(c(0,1), c(298.15, 298.15)))
  model$set_states("Ci", 400)
  model$set_states("Cc", 400)
  model$set_states("Ccyt", 400)
  model$set_states("PR", 100)
  model$set_states("gsw", 0.12)
  
  # Assign timepoints for the simulation
  model$set_time(c(1,3600))
  
  # Generate the simulation and watch for errors in the integration
  initial = try(cvode(model)[,names(model$States$Values)], silent = F)
  
  # Manage error
  if(inherits(initial,"try-error")) {print("error in model eval"); return(NULL)}
  
  # 2. Run simulation that corresponds to the experiment
  
  # Assign values of states variables achieved at the end of acclimation
  model$set_states(names(model$States$Values), initial[2,])
  
  # Assign the forcings corresponding to this measurement
  deltat = 15*60
  nt = 30
  irradiance = sort(rep(seq(0,2000,length.out = nt), 2))
  timepoints = (1:nt)*deltat
  timepoints = c(0,sort(c(timepoints, timepoints[-nt] + 1)))
  
  model$set_forcings("Ib", cbind(timepoints, irradiance*0.9))
  model$set_forcings("Ir", cbind(timepoints, irradiance*0.1))
  
  # Assign timepoints for the simulation
  model$set_time(0:max(timepoints))
  
  # Generate the simulation and watch for errors in the integration
  simulation = try(cvode(model), silent = F)
  
  # If there is an error, don't bother to do the rest of the simulation
  if(inherits(simulation,"try-error")) {print("error in model eval"); return(NULL)}
  
  # If everything went well, return list with photosynthesis and gs
  simulation[(1:nt)*deltat - 1,]
}

# Create the model object
model = generate_MiniModel_model()
# Setting for simulation
model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                   c(1e-10,1e-6,1e4,20,20,20, -1e-6))
model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
# Assign parameters from literature to the model
filter = which(names(parameters) %in% names(model$Parameters$Values))
model$set_parameters(names(parameters[filter]),unname(parameters[filter]))

# Run simulation
LRC = run_lrc(model)

rel = function(x) (x - min(x))/(max(x) - min(x))

# Plot values of fR, fRB and gs relative to maximum
png("Output/LRC.png", width = 8, height = 6, pointsize = 8, units = "cm", 
    res = 600, bg = "white", antialias = "default")
par(mfrow = c(1,1), xaxs = "i", yaxs = "i", las = 1, mar= c(4.0,4.0,0.5,1), mgp = c(2.3,1,0))
plot(c(0,2000), c(0,1), type = "n", xlab = expression(I~(mu*mol~m^{-2}~s^{-1})),
     ylab = "Relative value", ylim = c(0,1.05))
lines(LRC[,"PAR"], rel(LRC[,"Photo"]))
lines(LRC[,"PAR"], rel(LRC[,"fRB"]), col = 2, lty = 2)
lines(LRC[,"PAR"], rel(LRC[,"fR"]), col = 3, lty = 3)
lines(LRC[,"PAR"], rel(LRC[,"Cond"]), col = 4, lty = 4)
legend("bottomright", c("A", "RB", "R", "gs"), col = 1:4, lty = 1:4, bty = "n")
dev.off()







#