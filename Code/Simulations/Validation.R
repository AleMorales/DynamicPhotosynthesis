
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



# Square waves ------------------------------------------------------------

# Load the data
load("Intermediate/lightflecks.RData")
inputs = subset(lf_data_mean, Genotype %in% c("col", "npq1", "npq4", "rca2", "rwt43", "spsa") & Obs <= 1262) %>% arrange(Genotype, Amplitude, Obs)
lf = filter(inputs, Obs > 2) %>% mutate(Time = (Time - 3600)/60)
rm(lf_data_mean, lf_data_se)

# Functions to simulate the square lightflecks
run = function(model, data) {
  
  # 1. Run adaptation during acclimation
  
  # Assign constant conditions during acclimation
  model$set_forcings("Ib", cbind(data[1:2,"Time"],data[1:2,"PAR"]*0.1))
  model$set_forcings("Ir", cbind(data[1:2,"Time"],data[1:2,"PAR"]*0.9))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("H2OR", data[1:2,c("Time","H2OR")])
  model$set_forcings("CO2R", data[1:2,c("Time","CO2R")])
  model$set_forcings("Ta", data[1:2,c("Time","Tair")])
  model$set_forcings("Tl", data[1:2,c("Time","Tleaf")])
  model$set_states("Ci", data[1,"CO2R"])
  model$set_states("Cc", data[1,"CO2R"])
  model$set_states("Ccyt", data[1,"CO2R"])
  model$set_states("PR", 100)
  
  # Assign timepoints for the simulation
  model$set_time(c(1,7200))
  
  # Generate the simulation and watch for errors in the integration
  initial = try(cvode(model)[,names(model$States$Values)], silent = F)
  
  # Manage error
  if(inherits(initial,"try-error")) {print("error in model eval"); return(NULL)}
  
  # 2. Run simulation that corresponds to the experiment
  
  # Assign values of states variables achieved at the end of acclimation
  model$set_states(names(model$States$Values), initial[2,])
  
  # Assign the forcings corresponding to this measurement
  model$set_forcings("Ib", cbind(data[-(1:2),"Time"],data[-(1:2),"PAR"]*0.1))
  model$set_forcings("Ir", cbind(data[-(1:2),"Time"],data[-(1:2),"PAR"]*0.9))
  model$set_forcings("H2OR", data[-(1:2),c("Time","H2OR")])
  model$set_forcings("CO2R", data[-(1:2),c("Time","CO2R")])
  model$set_forcings("Ta",  data[-(1:2),c("Time","Tair")])
  model$set_forcings("Tl",  data[-(1:2),c("Time","Tleaf")])
  
  # Assign timepoints for the simulation
  model$set_time(data[-(1:2),"Time"])
  
  # Generate the simulation and watch for errors in the integration
  simulation = try(cvode(model), silent = F)
  
  # If there is an error, don't bother to do the rest of the simulation
  if(inherits(simulation,"try-error")) {print("error in model eval"); return(NULL)}
  
  # If everything went well, return list with photosynthesis and gs
  simulation[1:1260,c("A", "Photo")]
}
run_all = function() {
  
  listOutput = vector("list", 6)
  names(listOutput) = c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")
  
  for(G in names(listOutput)) {
    listOutput[[G]] = vector("list", 3)
    k = 1
    for(amp in c(100, 200, 500)) {
      # Create the basic model object and populate with parameters depending on the genotype
      if(G == "rwt43") {
        model = generate_MiniModelrwt43_model()
      } else {
        model = generate_MiniModel_model()
      }
      
      # Setting for simulation
      model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                         c(1e-10,1e-6,1e4,20,20,20, -1e-6))
      model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
      
      # Assign parameters from literature to the model
      filter = which(names(parameters) %in% names(model$Parameters$Values))
      model$set_parameters(names(parameters[filter]),unname(parameters[filter]))

      # Assign parameters specific to genotypes
      switch(G, 
             rca2 = {model$set_parameters(c("RCA", "Jmax25"), parameters_mutants[c("RCA_rca2", "Jmax25_rca2")])},
             spsa = {model$set_parameters("TPU25", parameters_mutants[["TPU25_spsa"]])},
             npq4 = {model$set_parameters(c("Jmax25", "KiqEp"), c(parameters_mutants["Jmax25_npq4"], 0)); 
               model$set_states("fP", 0)},
             npq1 = {model$set_parameters("KiqEz", 0); model$set_states("fZ", 0)},
             rwt43 = {model$set_parameters("Jmax25", parameters_mutants["Jmax25_rwt43"])}
      )
      listOutput[[G]][[k]] = run(model, subset(inputs, Genotype == G & Amplitude == amp))
      k = k + 1
    }
  }
  
  # Extract the individual values grouped by genotype and then amplitude
  output = matrix(NA, ncol = ncol(listOutput[[1]][[1]]), nrow = nrow(lf))
  pos = 1
  for(i in 1:length(listOutput)) {
    for(j in 1:length(listOutput[[i]])) {
      data = listOutput[[i]][[j]]
      newpos = pos + nrow(data) - 1
      output[pos:newpos,] = data
      pos = newpos + 1
    }
  }
  output
}

# Simulate all lightflecks and add simulated values
simLF = run_all()
lf %<>% group_by() %>% mutate(simPhoto = simLF[,2], simAn = simLF[,1],
                              Period = rep(rep(c(-1,120,60), times = c(60,600,600)), 18),
                              Direction = rep(rep(c(NA, rep(c("up", "down"),5), rep(c("up", "down"),10)),
                                                      times = c(60,rep(60,10), rep(30,20))), 18))

# Avg measured and simulated Photo per genotype, amplitude and period, including
lfAvg = group_by(lf, Genotype, Amplitude, Period) %>% summarise(Photo = mean(Photo), simPhoto = mean(simPhoto))
lf0 = filter(lfAvg, Period == -1) %>% group_by()
save(lf0, file = "Intermediate/lf0.RData")

# loess-based predictions of steady-state photosynthesis from LRC
load("Intermediate/LRC.RData")
LRC = rbind(LRC, mutate(lf0, PAR = 300) %>% select(-Amplitude, -Period))
modPhoto = plyr::dlply(LRC, "Genotype", function(x) loess(Photo~PAR, data = x))
modsimPhoto = plyr::dlply(LRC, "Genotype", function(x) loess(simPhoto~PAR, data = x))
LRC = group_by(LRC, Genotype) %>% mutate(fitPhoto = predict(modPhoto[[Genotype[1]]], newdata = data.frame(PAR = PAR)),
                                         fitsimPhoto = predict(modsimPhoto[[Genotype[1]]], newdata = data.frame(PAR = PAR)))
LRC = mutate(LRC, resid = Photo - fitPhoto, simresid = simPhoto - fitsimPhoto)
LRC %>% group_by(Genotype) %>% summarise(nRMSE = sd(resid)/mean(Photo), simRMSE = sd(simresid)/mean(simPhoto))
lf = group_by(lf, Genotype) %>% mutate(Photoss = predict(modPhoto[[Genotype[1]]], newdata = data.frame(PAR = PAR)),
                                       simPhotoss = predict(modsimPhoto[[Genotype[1]]], newdata = data.frame(PAR = PAR)))

# Calculate LFUE for whole cycles, and up/down cycles
LFUE = lfAvg = group_by(lf, Genotype, Amplitude, Period) %>% filter(Period > 0) %>% 
  summarise(LFUE = mean(Photo)/mean(Photoss), simLFUE = mean(simPhoto)/mean(simPhotoss),
            simLFUEAn = mean(simAn)/mean(simPhotoss))

LFUEcycle = lfAvg = group_by(lf, Genotype, Amplitude, Period, Direction) %>% 
  filter(Period > 0) %>% summarise(LFUE = mean(Photo)/mean(Photoss), simLFUE = mean(simPhoto)/mean(simPhotoss))

# Effect of delays on LFUE calculations
with(LFUE, quantile((simLFUEAn - simLFUE)/simLFUEAn*100, c(0.05,0.95)))
with(lf, quantile((simAn - simPhoto)/simAn*100, c(0.05,0.95)))


statsTotal = LFUEcycle %>% group_by(Genotype) %>%
    summarise(RMSE = sd(LFUE - simLFUE), R2 = 100*(1 - var(LFUE - simLFUE)/var(LFUE)))

statsTotalUp = LFUEcycle %>% filter(Direction == "up") %>% group_by(Genotype) %>%
  summarise(RMSE = sd(LFUE - simLFUE), R2 = 100*(1 - var(LFUE - simLFUE)/var(LFUE)))

statsTotalDow = LFUEcycle %>% filter(Direction == "down") %>% group_by(Genotype) %>%
  summarise(RMSE = sd(LFUE - simLFUE), R2 = 100*(1 - var(LFUE - simLFUE)/var(LFUE)))

# Validation plot -------------------------------------------------------------------------------------------------

png("Output/figureMeasuredModelledLFUE.png", width = 6, height = 18, pointsize = 11, units = "cm",
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(3,1), yaxs = "i", xaxs = "i", mar = c(4.0,4.5,0.5,0.8), las = 1, mgp = c(2.7,1,0), yaxs = "i", xaxs = "i")

# Full LFUE
plot(1,1, xlim = c(0.8,1.05), ylim = c(0.8,1.05), type = "n", xlab = expression(Measured~LFUE),
     ylab = expression(Simulated~LFUE))
abline(a = 0, b = 1)
genotype = c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")
for(i in 1:6)
  with(filter(LFUE, Genotype == genotype[i]), points(LFUE, simLFUE, col = i, pch = rep(1:3, each = 2)))
text(0.81,1.04,"A")
pos = legend("top", c(expression("Col-0 (100, 200, 500)"),
                      expression(italic("npq1-2")*" (100, 200, 500)"),
                      expression(italic("npq4-1")*" (100, 200, 500)"),
                      expression(italic("rwt43")*" (100, 200, 500)"),
                      expression(italic("rca-2")*" (100, 200, 500)"),
                      expression(italic("spsa1")*" (100, 200, 500)")), bty = "n", xjust = 0, ncol = 1, cex = 0.7)
points(x = rep(pos$text$x, each = 3)*c(0.97,0.98,0.99), 
       y = rep(pos$text$y, each = 3), 
       pch = rep(1:3, times = 6), 
       col = rep(1:6, each = 3))

# LFUE up-cycles
plot(1,1, xlim = c(0.5,1.1), ylim = c(0.5,1.1), type = "n", xlab = expression(Measured~"half-cycle"~LFUE),
     ylab = expression(Simulated~"half-cycle"~LFUE))
abline(a = 0, b = 1)
for(i in 1:6)
  with(filter(LFUEcycle, Direction == "up", Genotype == genotype[i]), points(LFUE, simLFUE, col = i, pch = rep(1:3, each = 2)))
text(0.52,1.08,"B")

# LFUE down-cycles
plot(1,1, xlim = c(1,2.7), ylim = c(1,2.7), type = "n", xlab = expression(Measured~"half-cycle"~LFUE),
     ylab = expression(Simulated~"half-cycle"~LFUE))
abline(a = 0, b = 1)
for(i in 1:6)
  with(filter(LFUEcycle, Direction == "down", Genotype == genotype[i]), points(LFUE, simLFUE, col = i, pch = rep(1:3, each = 2)))
text(1.05,2.65,"C")

dev.off()
