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
parnames = row.names(parameters)
parameters = parameters[,"value"]
names(parameters) = parnames
parameters["gamma3"]  = 1 - parameters["gamma2"] - parameters["gamma1"]
rm(data_summary, params_experiment, params_literature, experiment_summary)

run = function(model, parameters, PAR, PARtime) {
  
  # CVODE settings
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                     c(1e-10,1e-6,1e4,20,20,20, -1e-6))
  model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  
  # Assign parameters from literature to the model
  filter = which(names(parameters) %in% names(model$Parameters$Values))
  model$set_parameters(names(parameters[filter]),unname(parameters[filter]))
  
  
  # Simulate a LiCOR with 10% red - 90% blue
  P1 = PAR[1]
  model$set_forcings("Ib", cbind(c(0,1), c(P1,P1)/3))
  model$set_forcings("Ir", cbind(c(0,1), c(P1,P1)/3))
  model$set_forcings("Ig", cbind(c(0,1), c(P1,P1)/3))
  model$set_forcings("CO2R", cbind(c(0,1), c(400,400)))
  model$set_forcings("H2OR", cbind(c(0,1), c(20,20)))
  model$set_forcings("Ta", cbind(c(0,1), c(298.15, 298.15)))
  model$set_forcings("Tl", cbind(c(0,1), c(298.15, 298.15)))
  model$set_states("Ci", 400)
  model$set_states("Cc", 400)
  model$set_states("Ccyt", 400)
  tryCatch(model$set_states("PR", 25), error = function(x) NULL)
  
  # Calculate steady-state
  model$set_time(c(0,3600))
  steadyState = cvode(model)[2,names(model$States$Values)]
  model$set_states(names(steadyState), steadyState)
  
  # Simulate Transient
  model$set_forcings("Ib", cbind(PARtime, PAR/3))
  model$set_forcings("Ir", cbind(PARtime, PAR/3))
  model$set_forcings("Ig", cbind(PARtime, PAR/3))
  model$set_time(seq(PARtime[1], PARtime[length(PARtime)], 1))
  cvode(model)[,"A"]
}

# Simulate cloudflecks sequence

# Increase leaf photosynthesis and stomatal conductance for top leaf
cloudparameters = parameters
cloudparameters["Jmax25"] = parameters["Jmax25"]*3
cloudparameters["RB"] = parameters["RB"]*3
cloudparameters["TPU25"] = parameters["TPU25"]*3
cloudparameters["Rm25"] = parameters["Rm25"]*3
cloudparameters["gswm"] = parameters["gswm"]*3


cloudPAR = read.csv("cloudflecks_curve.csv")

control = run(generate_MiniModel_model(), cloudparameters, cloudPAR$Qt*0.45*4.57, (1:650)*60)[(350*60):(650*60)]

png(file = "CloudfleckSequence.png", width = 4, height = 4, units = "in",
    bg = "transparent", res = 1000)
par(mar = c(0.5,0.5,0.5,0.5), bg = "transparent")
plot(control, t = "l", lwd = 2, col = "darkgreen", xaxt = "n",
                   yaxt = "n", bg = "transparent", bty = "n")
dev.off()


# Simulate sunflecks sequence
sunPAR = read.csv("wheat.csv")

control = run(generate_MiniModel_model(), parameters, sunPAR$Below, 1:22497)[550:1200]

png(file = "SunfleckSequence.png", width = 4, height = 4, units = "in",
    bg = "transparent", res = 1000)
par(mar = c(0.5,0.5,0.5,0.5), bg = "transparent")
plot(control, t = "l", lwd = 2, col = "darkgreen", xaxt = "n",
     yaxt = "n", bg = "transparent", bty = "n")
dev.off()
