
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



# Induction states, gs and NPQ for light transients -----------------------------------

# Retrieve the induction data
load("Intermediate/transients.RData")

# Calculate SI for each transient
transients = subset(transients_mean, Genotype %in% c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")  &
                      Obs > 2) %>% arrange(Genotype, TransientType, Obs)
transients_se = subset(transients_se, Genotype %in% c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")  &
                         Obs > 2) %>% arrange(Genotype, TransientType, Obs)
modPhoto = transients %>% group_by(Genotype, TransientType) %>%
  do(fitPhoto = loess(Photo~Time, data= ., span = 0.05))
modCond = transients %>% group_by(Genotype, TransientType) %>%
  do(fitCond = loess(Cond~Time, data= ., span = 0.05))
sPhoto = unlist(lapply(modPhoto$fitPhoto, function(x) predict(x)))
sCond = unlist(lapply(modCond$fitCond, function(x) predict(x)))
transients = transients %>% group_by() %>% mutate(sPhoto = sPhoto, sCond = sCond) %>% group_by() %>% 
  mutate(NPQse = transients_se[["NPQ"]], Condse = transients_se[["Cond"]],
         Photose = transients_se[["Photo"]])

# Small differences in measurement conditions no significant effect on output
# Thus, use average conditions and one simulation per transient type
inputs = subset(transients_mean, Genotype %in% c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")) %>% arrange(Genotype, TransientType, Obs)

rm(transients_all, transients_se, transients_mean, modPhoto, sPhoto)

# Function to simulate an induction curve given a model object and the input data
run_transient = function(model, data) {
  
  # 1. Run adaptation during acclimation
  
  # Assign constant conditions during acclimation
  model$set_forcings("Ib", cbind(data[1:2,"Time"],data[1:2,"PAR"]*0.1))
  model$set_forcings("Ir", cbind(data[1:2,"Time"],data[1:2,"PAR"]*0.9))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("H2OR", data[1:2,c("Time","H2OR")])
  model$set_forcings("CO2R", data[1:2,c("Time","CO2R")])
  model$set_forcings("Ta", data[1:2,c("Time","Tair")])#cbind(1:2,rep(296.4,2)))
  model$set_forcings("Tl", data[1:2,c("Time","Tleaf")])#cbind(1:2,rep(296.4,2)))
  model$set_states("Ci", data[1,"CO2R"])
  model$set_states("Cc", data[1,"CO2R"])
  model$set_states("Ccyt", data[1,"CO2R"])
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
  model$set_forcings("Ib", cbind(data[-(1:2),"Time"],data[-(1:2),"PAR"]*0.1))
  model$set_forcings("Ir", cbind(data[-(1:2),"Time"],data[-(1:2),"PAR"]*0.9))
  model$set_forcings("H2OR", data[-(1:2),c("Time","H2OR")])
  model$set_forcings("CO2R", data[-(1:2),c("Time","CO2R")])
  model$set_forcings("Ta", data[-(1:2),c("Time","Tair")])
  model$set_forcings("Tl", data[-(1:2),c("Time","Tleaf")])
  
  # Assign timepoints for the simulation
  model$set_time(data[-(1:2),"Time"])
  
  # Generate the simulation and watch for errors in the integration
  simulation = try(cvode(model), silent = F)
  
  # If there is an error, don't bother to do the rest of the simulation
  if(inherits(simulation,"try-error")) {print("error in model eval"); return(NULL)}
  
  # If everything went well, return list with photosynthesis and gs
  simulation[,c("Photo", "Cond","NPQ","qE","qM","qI")]
}

run_all = function() {
  
  listOutput = vector("list", 6)
  names(listOutput) = c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")
  
  for(G in names(listOutput)) {
    listOutput[[G]] = vector("list", 5)
    names(listOutput[[G]]) = paste0("transient", 1:5)
    for(t in paste0("transient", 1:5)) {
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
      
      induction = run_transient(model, subset(inputs, Genotype == G & TransientType == t))
      listOutput[[G]][[t]] = induction
    }
  }
  
  # Extract the individual values grouped by genotype and then transientype
  output = matrix(NA, ncol = ncol(listOutput[[1]][[1]]), nrow = nrow(transients))
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

# Simulate all transients and retrieve photo, npq and cond
simTransients = run_all()
transients = transients %>% group_by() %>% mutate(simPhoto = simTransients[,1], simCond = simTransients[,2],
                                                  simNPQ = simTransients[,3], qE = simTransients[,4],
                                                  qM = simTransients[,5], qI = simTransients[,6])

# Plot NPQ for Col-0 and rwt43
ggplot(filter(transients, Genotype %in% c("col","rwt43")), aes(Time, qE, col = Genotype)) + geom_line() + 
  facet_wrap(~TransientType, ncol = 2)

PAR0 = c(transient1 = 0, transient2 = 70, transient3 = 800, transient4 = 130, transient5 = 600)
PARf = c(transient1 = 1000, transient2 = 800, transient3 = 130, transient4 = 600, transient5 = 200)
LRC = transients %>% group_by(Genotype, TransientType)  %>%  
  summarise(Photo0 = mean(Photo[1:30]), simPhoto0 = mean(simPhoto[1:30]),
            Photof = mean(Photo[(length(Photo) - 30):(length(Photo))]),
            simPhotof = mean(simPhoto[(length(simPhoto) - 30):(length(simPhoto))]),
            PAR0 = PAR0[TransientType[1]], PARf = PARf[TransientType[1]]) %>% group_by()
LRC = rbind(select(LRC, -Photof, -simPhotof, -PARf), select(LRC, -Photo0, -simPhoto0, -PAR0) %>% 
              rename(Photo0 = Photof, simPhoto0 = simPhotof, PAR0 = PARf)) %>%
  rename(PAR = PAR0, Photo = Photo0, simPhoto = simPhoto0) %>% select(-TransientType)
save(LRC, file = "Intermediate/LRC.RData")

calcIS = function(x) {
  xA = mean(x[1:60])
  xB = mean(x[(length(x) - 60):length(x)])
  (x - min(xA, xB))/abs(xB - xA)
}

calcIS2 = function(x) {
  m = which.min(x)
  xmin = mean(x[(m - 5):(m + 5)])
  xf = mean(x[(length(x) - 60):length(x)])
  (x - xmin)/abs(xf - xmin)
}
# Calculate relative Photo and gs from measurements and model
transients = transients %>% group_by(Genotype, TransientType)  %>% 
  mutate(IS = ifelse(TransientType %in% paste0("transient", c(1,2, 4)), calcIS(sPhoto), calcIS2(sPhoto)), 
         simIS = ifelse(TransientType %in% paste0("transient", c(1,2, 4)), calcIS(simPhoto), calcIS2(simPhoto)),
         ISgs = calcIS(sCond), simISgs = calcIS(simCond))

# Calculate for every genotype and TransientType %in% c(1, 2, 4) the values of t30, t50, t70 and t90
tableTransientUp = filter(transients, TransientType %in% paste0("transient", c(1,2,4))) %>% 
  group_by(Genotype, TransientType) %>% summarise(t30 = (Time[which.min(abs(IS - 0.3))] - Time[1])/60,
                                                  simt30 = (Time[which.min(abs(simIS - 0.3))] - Time[1])/60,
                                                  t50 = (Time[which.min(abs(IS - 0.5))] - Time[1])/60,
                                                  simt50 = (Time[which.min(abs(simIS - 0.5))] - Time[1])/60,
                                                  t70 = (Time[which.min(abs(IS - 0.7))] - Time[1])/60,
                                                  simt70 = (Time[which.min(abs(simIS - 0.7))] - Time[1])/60,
                                                  t90 = (Time[which.min(abs(IS - 0.90))] - Time[1])/60,
                                                  simt90 = (Time[which.min(abs(simIS - 0.90))] - Time[1])/60) %>%
  arrange(Genotype, TransientType)

calcDelta = function(x) {
  m = which.min(x)
  xmin = mean(x[(m - 5):(m + 5)])
  xf = mean(x[(length(x) - 60):length(x)])
  abs(xf - xmin)
}

# Calculate for every genotype and TransientType %in% c(4,5) the values of tmin, t50, t70 and t90 (the last 3 for recovery)
tableTransientDown = filter(transients, TransientType %in% paste0("transient", c(3,5))) %>% 
  group_by(Genotype, TransientType) %>%
  summarise(tmin = (Time[which.min(IS)] - Time[1])/60,
            simtmin = (Time[which.min(simIS)] - Time[1])/60,
            t50 = (Time[which.min(abs(IS[which.min(IS):length(IS)] - 0.5)) + which.min(IS)] - Time[1])/60,
            simt50 = (Time[which.min(abs(simIS[which.min(simIS):length(simIS)] - 0.5)) + which.min(simIS)] - Time[1])/60,
            t70 = (Time[which.min(abs(IS[which.min(IS):length(IS)] - 0.7)) + which.min(IS)] - Time[1])/60,
            simt70 = (Time[which.min(abs(simIS[which.min(simIS):length(simIS)] - 0.7)) + which.min(simIS)] - Time[1])/60,
            t90 = (Time[which.min(abs(IS[which.min(IS):length(IS)] - 0.9)) + which.min(IS)] - Time[1])/60,
            simt90 = (Time[which.min(abs(simIS[which.min(simIS):length(simIS)] - 0.9)) + which.min(simIS)] - Time[1])/60,
            DeltaPhoto = calcDelta(sPhoto), DeltaSimPhoto = calcDelta(simPhoto), Deltares = DeltaPhoto - DeltaSimPhoto)

obs = c(with(tableTransientUp, c(t30, t50, t70, t90)),with(tableTransientDown, c(tmin, t50, t70, t90)))
residuals = c(with(tableTransientUp, c(simt30 - t30, simt50 - t50, simt70 - t70, simt90 - t90)),
              with(tableTransientDown, c(simtmin - tmin, simt50 - t50, simt70 - t70, simt90 - t90)))

# Calculate R2 and RMSE per type of up transient and genotype
tableTransientUp = tableTransientUp %>% mutate(rest30 = simt30 - t30, rest50 = simt50 - t50, 
                                               rest70 = simt70 - t70, rest90 = simt90 - t90)
StatsTransient = tableTransientUp %>% group_by(Genotype, TransientType) %>%
                    summarise(RMSE = sd(c(rest30, rest50, rest70, rest90)),
                              R2 = (1 - var(c(rest30, rest50, rest70, rest90))/
                                        var(c(t30, t50, t70, t90)))*100)

# Calculate R2 and RMSE per type of down transient and genotype
tableTransientDown = tableTransientDown %>% mutate(restmin = simtmin - tmin, rest50 = simt50 - t50, 
                                               rest70 = simt70 - t70, rest90 = simt90 - t90)
StatsTransient = tableTransientDown %>% group_by(Genotype) %>%
  summarise(RMSE = sd(c(restmin, rest50, rest70, rest90)),
            R2 = (1 - var(c(restmin, rest50, rest70, rest90))/
                    var(c(tmin, t50, t70, t90)))*100)


plot_mutant = function(id, label, yaxis = F, xaxis = F, firstY = TRUE, firstX = TRUE) {
  with(subset(tableTransientUp, TransientType == id), {
    plot(t30, simt30, pch = 1:6, t = "p", ylim = c(1,30), xlim = c(1,30), col = 1,xaxt = "n", yaxt = "n", log = "xy")
    points(t50, simt50, pch = 1:6, col = 2)
    points(t70, simt70, pch = 1:6, col = 3)
    points(t90, simt90, pch = 1:6, col = 4) 
    lines(c(1,30),c(1,30), lty = 1)
  })
  if(yaxis) if(firstY) axis(2, c(1,2, 3,5,10,20, 30)) else axis(2, c(1,2, 3,5,10,20))
  if(xaxis) if(firstX) axis(1, c(1, 2, 3,5,10,20,30)) else axis(1, c(2, 3,5,10,20, 30))
  text(0.9,27,labels = label, pos = 4)
  
  pos = legend("top", c(expression(t[30]~"(Col-0, "*italic("npq1-2")*", "*italic("npq4-1")*","),
                        expression("      "~italic("spsa1")*", "*italic("rwt43")*", "*italic("rca-2")*")"),
                        expression(t[50]~"(idem)"),
                        expression(t[70]~"(idem)"), 
                        expression(t[90]~"(idem)")), 
               bty = "n", xjust = 0, ncol = 1, cex = 0.6)
  points(x = rep(seq(1.4,2.65,l = 6), each = 4), 
         y = rep(rev(c(14,16.5,19.5,26)), times = 6), 
         pch = rep(1:6, each = 4), 
         col = rep(1:4, times = 6), cex= 0.6) 
}

png("Output/figureMeasuredModelledIS.png", width = 14, height = 14, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(2,2), oma = c(4,4,1,1), mar = rep(0,4), las = 1, xpd = T, yaxs = "i", xaxs = "i")

plot_mutant("transient1", "A", T, F, T, F)
text(5,1.2, expression(0 %->% 1000))

plot_mutant("transient2", "B", F, F, F, F)
text(5,1.2, expression(70 %->% 800))

plot_mutant("transient4", "C", T, T, F, T)
text(5,1.2, expression(130 %->% 600))

symbols = rep(1:6, each = 2)
with(tableTransientDown, {
  # Plot kinetics
  plot(tmin, simtmin, ylim = c(1,30), xlim = c(1,30), log = "xy", pch = symbols, yaxt = "n", xaxt = "n",
       ylab = "", xlab = "")
  points(t50, simt50, col = 2, pch = symbols)
  points(t70, simt70, col = 3, pch = symbols)
  points(t90, simt90, col = 4, pch = symbols)
  lines(c(1,30),c(1,30), lty = 1)
  axis(1, c(2, 3,5,10,20, 30))
  
  pos = legend("top", c(expression(t["min"]~"(Col-0, "*italic("npq1-2")*", "*italic("npq4-1")*","),
                        expression("      "~italic("spsa1")*", "*italic("rwt43")*", "*italic("rca-2")*")"),
                        expression(t[50]~"(idem)"),
                        expression(t[70]~"(idem)"), 
                        expression(t[90]~"(idem)")), 
               bty = "n", xjust = 0, ncol = 1, cex = 0.6)
  points(x = rep(seq(1.4,2.65,l = 6), each = 4), 
         y = rep(rev(c(14,16.5,19.5,26)), times = 6), 
         pch = rep(1:6, each = 4), 
         col = rep(1:4, times = 6), cex= 0.6) 
  
  text(0.9,27,labels = "D", pos = 4)
})
mtext("Measured time index (min)", side = 1, line = 2.5, outer = T, cex = 0.85)
mtext("Simulated time index (min)", side = 2, line = 2.5, outer = T, las = 3, cex = 0.85)
text(5,1.5, expression(800 %->% 130))
text(5,1.2, expression(600 %->% 200))

dev.off()



tableTransientGs = filter(transients, TransientType %in% paste0("transient", 1:4)) %>%  
  group_by(Genotype, TransientType) %>% summarise(t30 = (Time[which.min(abs(ISgs - 0.3))] - Time[1])/60,
                                                  simt30 = (Time[which.min(abs(simISgs - 0.3))] - Time[1])/60,
                                                  t50 = (Time[which.min(abs(ISgs - 0.5))] - Time[1])/60,
                                                  simt50 = (Time[which.min(abs(simISgs - 0.5))] - Time[1])/60,
                                                  t70 = (Time[which.min(abs(ISgs - 0.7))] - Time[1])/60,
                                                  simt70 = (Time[which.min(abs(simISgs - 0.7))] - Time[1])/60,
                                                  t90 = (Time[which.min(abs(ISgs - 0.9))] - Time[1])/60,
                                                  simt90 = (Time[which.min(abs(simISgs - 0.9))] - Time[1])/60)

obs = with(tableTransientGs, c(t30, t50, t70, t90))
mod = with(tableTransientGs, c(simt30, simt50, simt70, simt90))

StatsTransient = tableTransientGs %>% group_by(Genotype) %>%
  summarise(RMSE = sd(c(t30, t50, t70, t90) - c(simt30, simt50, simt70, simt90)),
            R2 = (1 - var(c(t30, t50, t70, t90) - c(simt30, simt50, simt70, simt90))/
                    var(c(t30, t50, t70, t90)))*100)

# There is a crazy t90 value of 0.2 minutes...
tableTransientGs = filter(tableTransientUpGs, t90 > 1)

png("Output/figureMeasuredModelledgs.png", width = 8, height = 8, units = "cm", pointsize = 10,
    res = 600, bg = "white", antialias = "default")

symbols = as.numeric(sapply(tableTransientGs$TransientType, function(x) substr(x, 10,10)))
# Plot kinetics
par(mfrow = c(1,1), mar = c(3.8,3.8,0.5,1), las = 1, xaxs = "i", yaxs = "i", mgp = c(2.5,1,0))
with(tableTransientGs, {
  plot(t30, simt30, ylim = c(0,60), xlim = c(0,60), pch = symbols,
       ylab = "Simulated time index (min)", xlab = "Measured time index (min)")
  points(t50, simt50, col = 2, pch = symbols)
  points(t70, simt70, col = 3, pch = symbols)
  points(t90, simt90, col = 4, pch = symbols)
  abline(a = 0, b = 1, lty = 2)
  pos = legend("top", c(expression(t[30]~"(0 - 1000, 70 - 800"),
                        expression("      800 - 130, 130 - 600)"),
                        expression(t[50]~"(idem)"),
                        expression(t[70]~"(idem)"), 
                        expression(t[90]~"(idem)")), 
               bty = "n", xjust = 0, ncol = 1, cex = 0.8)
  points(x = rep(pos$text$x[-2], each = 4) - c(10,7,4,2), 
         y = rep(pos$text$y[-2], each = 4), 
         pch = rep(1:4, times = 4), 
         col = rep(1:4, each = 4), cex = 0.8)
})

dev.off()





npq_plotter = function(n) {
  gen = c("col", "npq4", "npq1", "spsa", "rwt43", "rca2")[n]
  with(subset(transients, Genotype == gen & flashon == 1), {
    Time = c(0,(Time - 3720)/60)
    points(Time, c(0,NPQ), col = n, pch = n)
    lines(Time, c(0,simNPQ), col = n, lty = n)
    errbar(Time[20], NPQ[20], NPQ[20] + 1.96*NPQse[20], NPQ[20] - 1.96*NPQse[20], add = T, pch = n, col = n, errbar.col = n)
    print(gen)
    print(sd(NPQ - simNPQ))
    print(100*(1 - var(NPQ - simNPQ)/var(NPQ)))
  })
}

q_plotter = function(n, var) {
  gen = c("col", "npq4", "npq1", "spsa", "rwt43", "rca2")[n]
  x = with(subset(transients, Genotype == gen & flashon == 1), c(0,(Time - 3720)/60))
  y = c(0, subset(transients, Genotype == gen & flashon == 1)[[var]])
  lines(x,y, col = n, lty = n)
}

# Plot NPQ and its components
png("Output/figureMeasuredModelledNPQ.png", width = 12, height = 12, units = "cm", pointsize = 10,
    res = 600, bg = "white", antialias = "default")
par(mfrow = c(2,2), mar = rep(0,4), oma = c(3.5,3.8,1,3.8), yaxs = "i", xaxs = "i", las = 1)

# Measured and modelled NPQ
plot(0,0, xlim = c(0,60), ylim = c(0,3), type = "b", xaxt = "n")
for(i in 1:6) npq_plotter(i)

mtext("NPQ", side = 2, las = 3, line = 2.7)
legend("bottomright",c(expression("Col-0"),expression(italic("npq4-1")),expression(italic("npq1-2")), expression(italic("spsa1")),
                       expression(italic("rwt43")),expression(italic("rca-2"))),
       col = c(1:6), pch = 1:6, lty = 1:6, bty = "n", ncol = 2, x.intersp = 1, y.intersp = 1, cex = 0.9)
text(2,2.9,"A")

# Modelled qE
plot(0,0, xlim = c(0,60), ylim = c(0,2.5), type = "b", xaxt = "n", yaxt = "n")
for(i in 1:6) q_plotter(i, "qE")
mtext("qE", side = 4, las = 3, line = 2.7)
axis(4, seq(0,2.5,0.5))
text(2,2.4,"B")

# Modelled qM
plot(0,0, xlim = c(0,60), ylim = c(0,0.35), type = "b", yaxt = "n", xaxt = "n")
for(i in 1:6) q_plotter(i, "qM")
axis(2, seq(0,0.30,0.05))
axis(1, seq(0,50,10))
mtext("qM", side = 2, las = 3, line = 2.7)
mtext("Time (min)", side = 1, line = 2.5)
text(2,0.34,"C")

# Modelled qI
plot(0,0, xlim = c(0,60), ylim = c(0,0.25), type = "b", yaxt = "n")
for(i in 1:6) q_plotter(i, "qI")
mtext("Time (min)", side = 1, line = 2.5)
mtext("qI", side = 4, las = 3, line = 2.7)
axis(4, seq(0,0.20,0.05))
text(2,0.24,"D")

dev.off()


# CO2 response curve ------------------------------------------------------

# A-Ci response curves from Kaiser et al. (2016)
load("Intermediate/aci.RData")
aci = subset(aci_df_genotype, Genotype %in% c("col", "npq1", "npq4", "rca2", "rwt43", "spsa") & CO2 > 50)
Ci = subset(aci_inputs_means, Genotype %in% c("col", "npq1", "npq4", "rca2", "rwt43", "spsa") & CO2R > 50)[,c("Genotype", "Time", "Ci")]
rm(aci_df, aci_df_genotype, aci_inputs_se, aci_inputs_means)


# Function to simulate an ACi curve given a model object and the input data with Ci measurements
run_aci = function(model, Ci) {
  
  # Timepoints at which steady-photo is logged
  index = c(3780, 3960, 4140, 4320, 5400, 5640, 5880, 6120) - 3600
  index = index[c(4,3,2,1,5,6,7,8)]
  
  # 1. Run adaptation during acclimation
  
  # Assign the forcings corresponding to this measurement
  model$set_forcings("Ib", cbind(1:2,rep(100,2)))
  model$set_forcings("Ir", cbind(1:2,rep(900,2)))
  model$set_forcings("Ig", cbind(1:2,c(0,0)))
  model$set_forcings("Ci", Ci[1:2,-1])
  model$set_forcings("Tl", cbind(1:2,rep(296.4,2)))
  
  # Assign timepoints for the simulation
  model$set_time(1:3600)
  
  # Generate the simulation and watch for errors in the integration
  initial = try(cvode(model)[,names(model$States$Values)])
  # Manage error
  if(inherits(initial,"try-error")) {print("error in model eval"); return(NULL)}
  
  #2.  Run actual simulation
  model$set_states(names(model$States$Values), initial[3600,])
  
  # Assign the forcings corresponding to this measurement
  model$set_forcings("Ci", Ci[-(1:2),-1])
  
  # Assign timepoints for the simulation
  model$set_time(3600:6120)
  
  # Generate the simulation and watch for errors in the integration
  simulations = try(cvode(model))
  
  # Manage error
  if(inherits(simulations,"try-error")) {print("error in model eval"); return(NULL)}
  
  # If everything went well, return list with the simulations
  simulations[index, "A"]
  
}

# Function to simulate all CO2 response curves in the orifinal experiment
run_all = function() {
  
  listOutput = vector("list", 6)
  names(listOutput) = c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")
  
  for(G in names(listOutput)) {
    # Create the basic model object and populate with parameters depending on the genotype
    if(G == "rwt43") {
      model = generate_MiniModelrwt43Ci_model()
    } else {
      model = generate_MiniModelCi_model()
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
    
    listOutput[[G]] = run_aci(model, subset(Ci, Genotype == G))
  }
  
  unlist(listOutput)
  
}

# Calculate simulated ACi
aci = aci %>% mutate(simA = run_all())

# Calculate steady-state gs
tableSteadyGsi = transients %>% group_by(Genotype, TransientType) %>% 
  summarise(Cond = mean(Cond[1:30]), Condse = mean(Condse[1:30]), 
            PAR = mean(PAR[1]), simCond = mean(simCond[1:30]))
tableSteadyGsf = transients %>% group_by(Genotype, TransientType) %>% 
  summarise(Cond = mean(Cond[(length(Cond) - 30):length(Cond)]), 
            Condse = mean(Condse[(length(Condse) - 30):length(Condse)]),
            PAR = mean(PAR[length(PAR)]), simCond = mean(simCond[(length(simCond) - 30):length(simCond)]))
tableSteadyGs = rbind(tableSteadyGsi, tableSteadyGsf) %>% arrange(Genotype, PAR)

# Calculate steady-state A-PAR
tableSteadyAi = transients %>% group_by(Genotype, TransientType) %>% 
  summarise(Photo = mean(Photo[1:30]), Photose = mean(Photose[1:30]), 
            PAR = mean(PAR[1]), simPhoto = mean(simPhoto[1:30]))
tableSteadyAf = transients %>% group_by(Genotype, TransientType) %>% 
  summarise(Photo = mean(Photo[(length(Photo) - 30):length(Photo)]), 
            Photose = mean(Photose[(length(Photose) - 30):length(Photose)]),
            PAR = mean(PAR[length(PAR)]), simPhoto = mean(simPhoto[(length(simPhoto) - 30):length(simPhoto)]))
tableSteadyA = rbind(tableSteadyAi, tableSteadyAf) %>% arrange(Genotype, PAR)


load("Intermediate/lf0.RData")
lf0 = filter(lf0, Amplitude == 100)  %>% 
  mutate(Photose = 0, TransientType = "NA", PAR = 300)

tableSteadyA = rbind(tableSteadyA %>% group_by(), lf0[,names(tableSteadyA)]) %>% arrange(Genotype, PAR)


# Generate R2 and RMSE for ACi curve
aci = aci %>% mutate(resPhoto = simA - Photo)
StatsAci = aci %>% group_by(Genotype) %>% summarise(R2 = (1 - var(resPhoto)/var(Photo))*100,
                                                    RMSE= sd(resPhoto))

# Generate R2 and RMSE for LRC
filter = rep(c(1,2,3,5,6,7, 9,11), times = 6) + rep(c(0, 11,22,33,44,54), each = 8)
LRC = tableSteadyA[filter,]
LRC = LRC %>% mutate(resPhoto = simPhoto - Photo)
StatsLRC = LRC %>% group_by(Genotype) %>% summarise(R2 = (1 - var(resPhoto)/var(Photo))*100,
                                                    RMSE= sd(resPhoto))
# Generate R2 and RMSE for LRC (gs)
filter = rep(c(1,2,3,5,6,8, 10), times = 6) + rep(c(0, 10,20,30,40,50), each = 7)
LRCCond = tableSteadyGs[filter,]
LRCCond = LRCCond %>% mutate(resCond = simCond - Cond)
StatsLRCCond = LRCCond %>% group_by(Genotype) %>% summarise(R2  = (1 - var(resCond)/var(Cond))*100,
                                                            RMSE = sd(resCond))


png("Output/figureMeasuredModelledSS.png", width = 7, height = 18, pointsize = 12, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(3,1), yaxs = "i", xaxs = "i", mar = c(3.9,4.5,0.5,1.5), las = 1, mgp = c(2.5,1,0))

# ACi 
with(subset(aci, Genotype == "col") , {
  plot(Ci, Photo, ylim = c(-1,25), xlim = c(0,1500),pch = 1, col = 1,t = "p",
       ylab = expression(italic(A)~(μmol~m^{-2}~s^{-1})), xlab = expression(italic(C["i"])~(μmol~mol^{-1})), 
       xaxt = "n", yaxt = "n")
  lines(Ci, simA)
  errbar(Ci[c(4,6,8)], Photo[c(4,6,8)], Photo[c(4,6,8)] + 1.96*Photo_se[c(4,6,8)], Photo[c(4,6,8)] - 1.96*Photo_se[c(4,6,8)], 
         add = T, pch = 1)
})
axis(1, at = seq(0,1500,500), tcl = -0.5)
axis(2, at = seq(0,25,5), tcl = -0.5)

with(subset(aci, Genotype == "npq1") , {
  points(Ci, Photo, xlab = "Time (min)", ylab = "Photo", col = 2, pch = 2)
  lines(Ci, simA, col = 2, lty = 2)
  errbar(Ci[c(4,6,8)], Photo[c(4,6,8)], Photo[c(4,6,8)] + 1.96*Photo_se[c(4,6,8)], Photo[c(4,6,8)] - 1.96*Photo_se[c(4,6,8)], 
         add = T, col = 2, errbar.col = 2, pch = 2)
})

with(subset(aci, Genotype == "npq4") , {
  points(Ci, Photo, xlab = "Time (min)", ylab = "Photo", col = 3, pch = 3)
  lines(Ci, simA, col = 3, lty = 3)
  errbar(Ci[c(4,6,8)], Photo[c(4,6,8)], Photo[c(4,6,8)] + 1.96*Photo_se[c(4,6,8)], Photo[c(4,6,8)] - 1.96*Photo_se[c(4,6,8)], 
         add = T, col = 3, errbar.col = 3, pch = 3)
})

with(subset(aci, Genotype == "rwt43") , {
  points(Ci, Photo, xlab = "Time (min)", ylab = "Photo", col = 4, pch = 4)
  lines(Ci, simA, col = 4, lty = 4)
  errbar(Ci[c(4,6,8)], Photo[c(4,6,8)], Photo[c(4,6,8)] + 1.96*Photo_se[c(4,6,8)], Photo[c(4,6,8)] - 1.96*Photo_se[c(4,6,8)], 
         add = T, col = 4, errbar.col = 4, pch = 4)
})
with(subset(aci, Genotype == "rca2") , {
  points(Ci, Photo, xlab = "Time (min)", ylab = "Photo", col = 5, pch = 5)
  lines(Ci, simA, col = 5, lty = 5)
  errbar(Ci[c(4,6,8)], Photo[c(4,6,8)], Photo[c(4,6,8)] + 1.96*Photo_se[c(4,6,8)], Photo[c(4,6,8)] - 1.96*Photo_se[c(4,6,8)], 
         add = T, col = 5, errbar.col = 5, pch = 5)
})
with(subset(aci, Genotype == "spsa") , {
  points(Ci, Photo, xlab = "Time (min)", ylab = "Photo", col = 6, pch = 6)
  lines(Ci, simA, col = 6, lty = 6)
  errbar(Ci[c(4,6,8)], Photo[c(4,6,8)], Photo[c(4,6,8)] + 1.96*Photo_se[c(4,6,8)], Photo[c(4,6,8)] - 1.96*Photo_se[c(4,6,8)], 
         add = T, col = 6, errbar.col = 6, pch = 6)
})
legend("bottomright",c(expression("Col-0"),expression(italic("npq1-2")), expression(italic("npq4-1")),
                       expression(italic("rwt43")),expression(italic("rca-2")), expression(italic("spsa1"))), 
       col = c(1:6), pch = 1:6, lty = 1:6, bty = "n", ncol = 1)
text(45,24,"A")

# APAR
genotypes = c("col", "npq1", "npq4", "rwt43", "rca2", "spsa")
filter = rep(c(1,2,3,5,6,7, 9,11), times = 6) + rep(c(0, 11,22,33,44,54), each = 8)
A_plotter = function(i) {
  with(subset(tableSteadyA[filter,], Genotype == genotypes[i]), {
    points(PAR, Photo, col = i, pch = i)
    lines(PAR, simPhoto, col = i, lty = i)
    errbar(PAR[c(4,6)], Photo[c(4,6)], Photo[c(4,6)] + 1.96*Photose[c(4,6)], Photo[c(4,6)] - 1.96*Photose[c(4,6)], 
           add = T, col = i, errbar.col = i, pch = i)
  })
}

plot(1,1,t = "n", xlim = c(-10,1100), ylim = c(-1,15), ylab = expression(italic(A)~(mu*mol~m^{-2}~s^{-1})),
     xlab = expression(I~(mu*mol~m^{-2}~s^{-1})))
A_plotter(1)
A_plotter(2)
A_plotter(3)
A_plotter(4)
A_plotter(5)
A_plotter(6)
text(35,14.3,"B")


# gs-PAR
filter = rep(c(1,2,3,5,6,8, 10), times = 6) + rep(c(0, 10,20,30,40,50), each = 7)
par(mgp = c(2.7,1,0))

gs_plotter = function(i) {
  with(subset(tableSteadyGs[filter,], Genotype == genotypes[i]), {
    points(PAR, Cond, col = i, pch = i)
    lines(PAR, simCond, col = i, lty = i)
    errbar(PAR[c(4,6)], Cond[c(4,6)], Cond[c(4,6)] + 1.96*Condse[c(4,6)], Cond[c(4,6)] - 1.96*Condse[c(4,6)], 
           add = T, col = i, errbar.col = i, pch = i)
  })
}


plot(1,1,t = "n", xlim = c(-10,1100), ylim = c(0,0.30), ylab = expression(g[s]~(mol~m^{-2}~s^{-1})),
     xlab = expression(I~(mu*mol~m^{-2}~s^{-1})))
gs_plotter(1)
gs_plotter(2)
gs_plotter(3)
gs_plotter(4)
gs_plotter(5)
gs_plotter(6)
text(35,0.285,"C")

dev.off()



# Supplemental Figure: Transient Decrease Irradiance --------------------------------------------------------------

transient3 = filter(transients, TransientType  == "transient3")
transient5 = filter(transients, TransientType  == "transient5")

selectData = function(data, genotype) {
  data = subset(data, Genotype == genotype)
  data = data[seq(1,1200, 5),]
}


png("Output/figureMeasuredModelledDecrease.png", width = 12, height = 18, pointsize = 12, units = "cm", 
    res = 600, bg = "white", antialias = "default")


par(mfrow = c(3,2), oma = c(4,4,0.5,0.5), mar = rep(0,4), yaxs = "i", xaxs = "i", 
    las = 1, mgp = c(2.5,1,0))

with(selectData(transient3, "col"), {
  plot((Time - Time[1])/60, Photo, xaxt = "n", xlab = "", ylim = c(0,15), xlim = c(0,20), yaxt = "n", 
       ylab = "", col = rgb(0,0,0,0.5))
  lines((Time - Time[1])/60, simPhoto, lwd = 2)
  axis(2, seq(0,15,3))
  mtext(side = 2, text = expression(italic(A)~(mu*mol~m^{-2}~s^{-1})), las = 3, line = 2.3, cex = 0.8)
  legend("topright", c(expression(800 %->% 130), expression(600 %->% 200)), col = 1:2, lty = 1:2, pch = 1:2, bty = "n")
})
with(selectData(transient5, "col"), {
  points((Time - Time[1])/60, Photo, col = rgb(1,0,0,0.5), pch = 2)
  lines((Time - Time[1])/60, simPhoto, col = 2, lty = 2, lwd = 2)
})
text(1.5,14.5,"Col-0")

with(selectData(transient3, "spsa"), {
  plot((Time - Time[1])/60, Photo, xaxt = "n", xlab = "", yaxt = "n", ylab = "", ylim = c(0,15), xlim = c(0,20), col = rgb(0,0,0,0.5))
  lines((Time - Time[1])/60, simPhoto, lwd = 2)
})
with(selectData(transient5, "spsa"), {
  points((Time - Time[1])/60, Photo, col = rgb(1,0,0,0.5), pch = 2)
  lines((Time - Time[1])/60, simPhoto, col = 2, lty = 2, lwd = 2)
})
text(1.8,14.5, expression(italic(spsa1)))

with(selectData(transient3, "rca2"), {
  plot((Time - Time[1])/60, Photo, xaxt = "n", xlab = "", ylim = c(0,15), xlim = c(0,20), yaxt = "n", 
       ylab = "", col = rgb(0,0,0,0.5))
  lines((Time - Time[1])/60, simPhoto, lwd = 2)
  axis(2, seq(0,12,3))
  mtext(side = 2, text = expression(italic(A)~(mu*mol~m^{-2}~s^{-1})), las = 3, line = 2.3, cex = 0.8)
})
with(selectData(transient5, "rca2"), {
  points((Time - Time[1])/60, Photo, col = rgb(1,0,0,0.5), pch = 2)
  lines((Time - Time[1])/60, simPhoto, col = 2, lty = 2, lwd = 2)
})
text(1.5,14.5, expression(italic("rca2")))

with(selectData(transient3, "rwt43"), {
  plot((Time - Time[1])/60, Photo, xaxt = "n", xlab = "", yaxt = "n", ylab = "", ylim = c(0,15), xlim = c(0,20),
       col = rgb(0,0,0,0.5))
  lines((Time - Time[1])/60, simPhoto, lwd = 2)
})
with(selectData(transient5, "rwt43"), {
  points((Time - Time[1])/60, Photo, col = rgb(1,0,0,0.5), pch = 2)
  lines((Time - Time[1])/60, simPhoto, col = 2, lty = 2, lwd = 2)
})
text(1.7,14.5, expression(italic(rwt43)))

with(selectData(transient3, "npq1"), {
  plot((Time - Time[1])/60, Photo, ylim = c(0,15), xlim = c(0,20), xaxt = "n", yaxt = "n", 
       ylab = "", col = rgb(0,0,0,0.5))
  lines((Time - Time[1])/60, simPhoto, lwd = 2)
  axis(1, seq(0,15,5))
  axis(2, seq(0,12,3))
  mtext(side = 2, text = expression(italic(A)~(mu*mol~m^{-2}~s^{-1})), las = 3, line = 2.3, cex = 0.8)
  mtext(side = 1, text = "Time (min)", las = 1, line = 2.3, cex = 0.8)
})
with(selectData(transient5, "npq1"), {
  points((Time - Time[1])/60, Photo, col = rgb(1,0,0,0.5), pch = 2)
  lines((Time - Time[1])/60, simPhoto, col = 2, lty = 2, lwd = 2)
})
text(1.8,14.5, expression(italic("npq1-2")))

with(selectData(transient3, "npq4"), {
  plot((Time - Time[1])/60, Photo, ylab = "", yaxt = "n", ylim = c(0,15), xlim = c(0,20), col = rgb(0,0,0,0.5))
  lines((Time - Time[1])/60, simPhoto, lwd = 2)
  mtext(side = 1, text = "Time (min)", las = 1, line = 2.3, cex = 0.8)
})
with(selectData(transient5, "npq4"), {
  points((Time - Time[1])/60, Photo, col = rgb(1,0,0,0.5), pch = 2)
  lines((Time - Time[1])/60, simPhoto, col = 2, lty = 2, lwd = 2)
})
text(1.8,14.5, expression(italic("npq4-1")))

dev.off()
