# Load libraries --------------------------------------------------------------------------------------------------
library(MiniModel)
library(nlfit)
library(dplyr)
library(ggplot2)
library(doParallel)


# Fit dynamic stomatal conductance ----------------------------------------

load(file = "Intermediate/ParametersExperiment.RData")

# Retrieve the data
load("Intermediate/transients.RData")

# Extract induction transients from col-0 and aba-2 (rest of genotypes have the same the gs as col-0)
transients = transients_mean %>% 
  filter(Genotype == "col" & TransientType %in% paste0("transient", 1:5)) %>% 
  arrange(TransientType, Obs)

# Light intensities use at the beginning and end of each transient
PAR0 = transients %>% group_by(TransientType) %>% summarise(PAR = PAR[1])
PARf = transients %>% group_by(TransientType) %>% summarise(PAR = PAR[length(PAR)])
VPD0 = transients %>% group_by(TransientType) %>% summarise(VPD = VPD[1])
VPDf = transients %>% group_by(TransientType) %>% summarise(VPD = VPD[length(VPD)])
GS0 = transients %>% group_by(TransientType) %>% summarise(gs0 = Cond[1])

fits_Cond = foreach(t = paste0("transient", 1:5)) %do% {
    # Select right curve
    data = subset(transients, TransientType == t & Obs > 1)
    # Estimate initial values for steady-state gs
    gs0 = data[1,"Cond"]
    gsm = data[nrow(data),"Cond"]
    # Remove first 2 minutes of measurement due to effect of open gas exchange system
    fit = nlfit(Cond~gs0 + (gsm - gs0)*(1 - exp(-Obs*k)),
                start = c(gs0 = gs0, gsm = gsm, k = 1/15/60),
                data = data[-(1:60),], lower = rep(0,3))
  }

names(fits_Cond) = paste0("transient", 1:5)

# Variability of Kgs among replicates and transients
Kgsi = c(fits_Cond$transient1$par["k"], fits_Cond$transient2$par["k"], fits_Cond$transient4$par["k"])
Kgsd = c(fits_Cond$transient3$par["k"], fits_Cond$transient5$par["k"])

# Calculate confidence intervals of rate constant for each genotype
if(.Platform$OS.type == "windows") {
  cl <- makeCluster(8) # Hyperthreading does not offset overhead of sockets
} else {
  cl <- makeForkCluster(8)
}
registerDoParallel(cl)

profile_Kgs = llply(fits_Cond, function(x) profile(x, cl = cl, maxsteps = 20, lower_delta = 0.01))
l_ply(profile_Kgs, function(x) plot(x))
confint_Kgs = llply(profile_Kgs, function(x) confint(x))

stopCluster(cl)

# Extract steady-state parameters for each genotype
Condss = as.data.frame(cbind(laply(fits_Cond, function(x) x$par), PAR0 = PAR0$PAR, PARf = PARf$PAR, VPD0 = VPD0$VPD, VPDf = VPDf$VPD))
Condss = transform(Condss, gs0 = GS0$gs0)
Condss = with(Condss, data.frame(PAR = c(PAR0, PARf), gs = c(gs0,gsm), VPD = c(VPD0,VPDf)))
Condss = Condss %>% arrange(PAR)

ggplot(Condss, aes(PAR, gs)) + geom_point()

# Fit steady-state model

create_model = function(PAR, VPD) {
  function(pars) {
    with(as.list(pars), {
      fI_a = thetafI
      fI_b = -(1 + fI0 + alphafI*PAR)
      fI_c = fI0 + alphafI*PAR
      fI = (-fI_b - sqrt(fI_b^2 - 4*fI_a*fI_c))/(2*fI_a)
      fvpd = 1/(1 + VPD/D0)
      return(fI*fvpd*gswm)
    })
  }
}

# With the data available it is not possible to alpha, theta and D0 simultaneously as:
# 1 - The low light region is not very sensitive to theta and D0 (so any value of theta or D0 would do)
# 2 - The high light region depends on both theta and D0
# The result is very strong correlation between the three parameters. So I fix theta and D0 to "reasonable" values

# Fit Col-0
model = create_model(Condss$PAR, Condss$VPD)
start = c(fI0 = 0.37, alphafI = 1e-3, thetafI = 0.9, D0 = 1, gswm = 1)

fitCol = nlfit(model = model, start = start, obs = Condss$gs, 
               lower = c(0.1, 1e-4, 0, 0.1, 0.3),
               upper = c(1,1e-1,1.0,2,4), fixed = c("D0", "thetafI"), 
               options = list(maxiter = 5e3, xrtol = 0, fatol = 1e-10, xatol = 0))

# Calculate confidence intervals
profile = profile(fitCol, tolerance = 10)
confintCol = confint(profile)



# Store the parameters with their confidence intervals in the table with parameter information
# Given the similarities, we choose one transient for increasing and one for decreasing irradiance

params_experiment["Kgsi", ] = confint_Kgs$transient1["k",]
params_experiment["Kgsd", ] = confint_Kgs$transient3["k",]
params_experiment["fI0", ] = confintCol["fI0",]
params_experiment["alphafI", ] = confintCol["alphafI",]
params_experiment["thetafI", ] = c(0.9, NA, NA)
params_experiment["gswm", ] = confintCol["gswm",]

experiment_summary["Kgs",] = c(paste0(c("Kgsi", "Kgsd"), collapse = ", "), 
                                   NRMSE = round(rmse(fits_Cond$transient1)/mean(fits_Cond$transient1$obs), 3), 
                                   CD = round(cd(fits_Cond$transient1), 3), 
                               n = length(c(fits_Cond$transient1$obs, fits_Cond$col$transient4$obs)))

experiment_summary["gss",] = c(paste0(c("fI0", "alphafI", "thetafI", "gswm"), collapse = ", "), 
                               NRMSE = round(rmse(fitCol)/mean(fitCol$obs), 3), 
                               CD = round(cd(fitCol), 3), n = length(fitCol$obs))

save(params_experiment, experiment_summary, file = "Intermediate/ParametersExperiment.RData")
rm(list = ls())



# Calculate dark respiration ----------------------------------------------

load(file = "Intermediate/ParametersExperiment.RData")

# Extract data from transients
load("Intermediate/transients.RData")
transients_all = subset(transients_all, Genotype == "col")

# Calculate respiration per genotype
Rdk = transients_all %>%
  filter(Obs %in% 1:120, TransientType == "transient1") %>%
  mutate(Resp = -Photo) %>% group_by(ID) %>% summarise(Resp = mean(Resp)) %>% group_by()


# There does not seem to be any significant differences across genotypes so we calculate a global value
Rdk23  = Rdk %>% summarise(value = mean(Resp), Rdk_se = sd(Resp)/sqrt(length(Resp))) %>%
  mutate(lower = value - 1.96*Rdk_se, upper = value + 1.96*Rdk_se) %>% select(value, lower, upper)

# Correct for temperature
Tleaf = transients_all %>%
  filter(Obs %in% 1:120, TransientType == "transient1") %>% summarise(Tleaf = mean(Tleaf))
Tleaf = Tleaf[[1]]
R = 8.31
Tref = 298.15
load("Intermediate/ParametersLiterature.RData")
RdHa = params_literature["DHaRm","value"]
Rd25 = Rdk23$value/exp((RdHa*1e3*(Tleaf - Tref))/(Tref*R*Tleaf))
f = Rd25/Rdk23$value
Rdk25 = Rdk23*f
row.names(Rdk25) = "Rm25"

params_experiment["Rm25", ] = Rdk25


# Construct the residuals
Res = Rdk$Resp - Rdk23$value
experiment_summary["Rm25",] = c("Rm25",
                                NRMSE = round(sd(Res)/mean(Rdk$Resp), 3),
                                CD = round(1 - sum(Res^2)/sum(Rdk$Resp^2), 3),
                                n = nrow(Rdk))

save(params_experiment, experiment_summary, file = "Intermediate/ParametersExperiment.RData")
rm(list = ls())

# Fit ACi and Induction curve ---------------------------------------------

load(file = "Intermediate/ParametersExperiment.RData")

# Retrieve the ACi data
# The simulations sometimes become numerically unstable at CO2R = 50 ppm, so we avoid the lowest point (still get good estimate of Vcmax)
load("Intermediate/aci.RData")
aci = subset(aci_df_genotype, Genotype %in% c("col", "npq1", "npq4", "rca2", "spsa", "rwt43") & CO2 > 50) %>% arrange(Genotype, Ci)
Ci = subset(aci_inputs_means, Genotype %in% c("col", "npq1", "npq4", "rca2", "spsa", "rwt43") & CO2R > 50)[,c("Genotype", "Time", "Ci")]
rm(aci_df, aci_df_genotype, aci_inputs_se, aci_inputs_means)

# Retrieve the induction data
load("Intermediate/transients.RData")

# Extract induction transients from col-0
# Fit loess to each transient in order to calculate relative photosynthetic rates during induction
# The reason is that there is variation in Vcmax across replicates, which affects absolute values
transients = subset(transients_mean, Genotype %in% c("col", "npq1", "npq4", "rca2", "spsa", "rwt43") & 
                      Obs > 2) %>% arrange(Genotype, TransientType, Obs) 
modPhoto = transients %>% group_by(Genotype, TransientType) %>%
  do(fitPhoto = loess(Photo~Time, data= ., span = 0.05))
sPhoto = unlist(lapply(modPhoto$fitPhoto, function(x) predict(x)))

modCond = transients %>% group_by(Genotype, TransientType) %>%
  do(fitCond = loess(Cond~Time, data= ., span = 0.05))
sCond = unlist(lapply(modCond$fitCond, function(x) predict(x)))

transients = transients %>% mutate(sPhoto = sPhoto, sCond = sCond) %>% group_by(Genotype, TransientType)  %>% 
  mutate(rPhoto = {photoA = mean(sPhoto[1:60])
  photoB = mean(sPhoto[(length(sPhoto) - 60):length(sPhoto)])
  (sPhoto - min(photoA, photoB))/abs(photoB - photoA)},
  rCond = {CondA = mean(sCond[1:60])
  CondB = mean(sCond[(length(sCond) - 60):length(sCond)])
  (sCond - min(CondA, CondB))/abs(CondB - CondA)})

inputs = subset(transients_mean, Genotype %in% c("col", "npq1", "npq4", "rca2", "spsa", "rwt43")) %>% arrange(Genotype, TransientType, Obs)
modCi = inputs %>% group_by(Genotype, TransientType) %>% do(fitCi = loess(Ci~Time, data= ., span = 0.05))
sCi = unlist(lapply(modCi$fitCi, function(x) predict(x)))
inputs = inputs %>% mutate(sCi = sCi)

rm(transients_all, transients_se, transients_mean, modPhoto, sPhoto, modCond, sCond)


# Function to simulate an ACi curve given a model object and the input data
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

# Function to simulate an induction curve given a model object and the input data
run_transient = function(model, data) {
  
  # 1. Run adaptation during acclimation
  
  # Assign constant conditions during acclimation
  model$set_forcings("Ib", cbind(data[1:2,"Time"],data[1:2,"PAR"]*0.1))
  model$set_forcings("Ir", cbind(data[1:2,"Time"],data[1:2,"PAR"]*0.9))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("H2OR", data[1:2,c("Time","H2OR")])
  model$set_forcings("CO2R", data[1:2,c("Time","CO2R")])
  model$set_forcings("Ta", cbind(1:2,rep(296.4,2)))
  model$set_forcings("Tl", cbind(1:2,rep(296.4,2)))
  model$set_states("Ci", data[1,"CO2R"])
  model$set_states("Cc", data[1,"CO2R"])
  model$set_states("Ccyt", data[1,"CO2R"])
  model$set_states("PR", 100)
  # Assign timepoints for the simulation
  model$set_time(c(1,45*60))
  
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
  
  # Assign timepoints for the simulation
  model$set_time(data[-(1:2),"Time"])
  
  # Generate the simulation and watch for errors in the integration
  simulation = try(cvode(model), silent = F)
  
  # If there is an error, don't bother to do the rest of the simulation
  if(inherits(simulation,"try-error")) {print("error in model eval"); return(NULL)}
  
  # If everything went well, return list with relative photosynthesis and NPQ
  Photo = simulation[,"Photo"]
  photoA = mean(Photo[1:60])
  photoB = mean(Photo[(length(Photo) - 60):length(Photo)])
  IS = (Photo - min(photoA, photoB))/abs(photoB - photoA)
      
  cbind(IS, simulation[,"NPQ"])
}


# Function to run all simulations and return the output of the simulations
# All predictions (An, IS and qE) are returned in a single vector
run_all = function(pars) {
  
  # Create the basic model object and populate with parameters depending on the genotype
  model = generate_MiniModel_model()
  modelCi = generate_MiniModelCi_model()
  model43 = generate_MiniModelrwt43_model()
  modelCi43 = generate_MiniModelrwt43Ci_model()
  
  # Setting for simulation
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                     c(1e-10,1e-10,1e4,20,20,20, -1e-6))
  model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  modelCi$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                       c(1e-10,1e-10,1e4,20,20,20, -1e-6))
  modelCi$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  model43$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                       c(1e-10,1e-10,1e4,20,20,20, -1e-6))
  model43$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  modelCi43$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                       c(1e-10,1e-10,1e4,20,20,20, -1e-6))
  modelCi43$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  
  # Assign parameters from literature to the model
  load("Intermediate/ParametersLiterature.RData")
  filter = which(!is.na(params_literature$value) & row.names(params_literature) %in% names(model$Parameters$Values))
  model$set_parameters(row.names(params_literature)[filter], params_literature[filter, "value"])
  
  load("Intermediate/ParametersExperiment.RData")
  filter = which(!is.na(params_experiment$value) & row.names(params_experiment) %in% names(model$Parameters$Values))
  model$set_parameters(row.names(params_experiment)[filter], params_experiment[filter, "value"])

  # Assign pars of wildtype
  gamma1 = pars[12]
  gamma2 = pars[13]
  if(gamma1 + gamma2 > 1) return(1e6)
  gamma3 = 1 - gamma1 - gamma2
  model$set_parameters(parnames, pars[parnames])
  model$set_parameters("gamma3", gamma3)
  

  # Assign parameters to the other models
  filter = which(names(model$Parameters$Values) %in% names(modelCi$Parameters$Values))
  parnames = names(model$Parameters$Values)[filter]
  modelCi$set_parameters(parnames, model$get_parameters(parnames))
  
  filter = which(names(model$Parameters$Values) %in% names(model43$Parameters$Values))
  parnames = names(model$Parameters$Values)[filter]
  model43$set_parameters(parnames, model$get_parameters(parnames))
  
  filter = which(names(model$Parameters$Values) %in% names(modelCi43$Parameters$Values))
  parnames = names(model$Parameters$Values)[filter]
  modelCi43$set_parameters(parnames, model$get_parameters(parnames))
  
  listOutput = try(foreach(G = c("col", "npq1", "npq4", "rca2", "rwt43", "spsa")) %dopar% {
    
    if(G == "rwt43") {
      local_model = model43$clone()
      local_modelCi = modelCi43$clone()
    } else {
      local_model = model$clone()
      local_modelCi = modelCi$clone()
    }


    # Assign parameters specific to genotypes
    switch(G, 
           rca2 = {local_model$set_parameters(c("RCA", "Jmax25"), pars[c("RCA_rca2", "Jmax25_rca2")])},
           spsa = {local_model$set_parameters("TPU25", pars[["TPU25_spsa"]])},
           npq4 = {local_model$set_parameters(c("Jmax25", "KiqEp"), c(pars["Jmax25_npq4"], 0)); 
             local_model$set_states("fP", 0)},
           npq1 = {local_model$set_parameters("KiqEz", 0); local_model$set_states("fZ", 0)},
           rwt43 = {local_model$set_parameters("Jmax25", pars["Jmax25_rwt43"])}
    )
    switch(G, 
           rca2 = {local_modelCi$set_parameters(c("RCA", "Jmax25"), pars[c("RCA_rca2", "Jmax25_rca2")])},
           spsa = {local_modelCi$set_parameters("TPU25", pars[["TPU25_spsa"]])},
           npq4 = {local_modelCi$set_parameters(c("Jmax25", "KiqEp"), c(pars["Jmax25_npq4"], 0)); 
             local_modelCi$set_states("fP", 0)},
           npq1 = {local_modelCi$set_parameters("KiqEz", 0); local_modelCi$set_states("fZ", 0)},
           rwt43 = {local_modelCi$set_parameters("Jmax25", pars["Jmax25_rwt43"])}
    )
    
    aci = run_aci(local_modelCi, subset(Ci, Genotype == G))
    if(is.null(aci)) stop("Error in aci for ", G)
    induction1 = run_transient(local_model, subset(inputs, Genotype == G & TransientType == "transient1"))
    if(is.null(induction1)) stop("Error in induction1 for ", G)
    induction2 = run_transient(local_model, subset(inputs, Genotype == G & TransientType == "transient2"))
    if(is.null(induction2)) stop("Error in induction2 for ", G)
    induction3 = run_transient(local_model, subset(inputs, Genotype == G & TransientType == "transient3"))
    if(is.null(induction3)) stop("Error in induction3 for ", G)
    induction4 = run_transient(local_model, subset(inputs, Genotype == G & TransientType == "transient4"))
    if(is.null(induction4)) stop("Error in induction4 for ", G)
    induction5 = run_transient(local_model, subset(inputs, Genotype == G & TransientType == "transient5"))
    if(is.null(induction5)) stop("Error in induction5 for ", G)
    
    # Make sure that we do not have a large slope got the last NPQ part of wildtype
    if(G == "col") {
      last = (nrow(induction1) - 300):nrow(induction1)
      time = seq_along(last)
      fit = lm(induction1[last,2]~time)
      NPQslope = coef(fit)[2]
      if(NPQslope > 0.5e-4) stop("Error Kinh for ", G)
    }
    
    list(aci = aci, induction = rbind(induction1, induction2, induction3, induction4, induction5))
  }, silent = TRUE)
  
  if(class(listOutput) != "list") return(rep(1e6, Nobs))

  # Extract the individual values
  aci = unlist(lapply(listOutput, function(x) x$aci))
  Photo = unlist(lapply(listOutput, function(x) x$induction[,1]))
  NPQ = unlist(lapply(listOutput, function(x) x$induction[,2]))
  filter = which(transients$NPQ > 0)
  NPQ = NPQ[filter]
  
  # Return vector of predictions
  c(aci, Photo, NPQ)
  
}



# Parameters and boundaries
pars = c(RB = 14.25, Jmax25 = 122, TPU25 = 7, PhiqEmax = 0.2, Kinh0 = 1e-7,fprot = 5e-7, 
         KiqEp = 0.116, KiqEz = 2.6e-3, Krca = 4.10e-05, fRBmin = 2.42e-1 , Vrmax = 60, 
         gamma1 = 0.4, gamma2 = 0.5, RCA_rca2 = 24.42, Jmax25_rca2 = 100, Jmax25_npq4 = 100, TPU25_spsa = 5.998177,
         Jmax25_rwt43 = 139.5547)

upper = c(RB = 21, Jmax25 = 220, TPU25 = 9, PhiqEmax = 1, Kinh0 = 1e-5,fprot = 1e-5, 
          KiqEp = Inf, KiqEz = 10,  Krca = 1e-3, fRBmin = 1, Vrmax = 200, gamma1 = 1, gamma2 = 1, RCA_rca2 = 100, 
          Jmax25_rca2 = 200, Jmax25_npq4 = 200, TPU25_spsa = 8, Jmax25_rwt43 = 150)

lower = c(RB = 10, Jmax25 = 100, TPU25 = 6, PhiqEmax = 0, Kinh0 = 1e-9,fprot = 1e-9,
          KiqEp = 0, KiqEz = 0,  Krca = 1e-6, fRBmin = 1e-10, Vrmax = 0, gamma1 = 0.1, gamma2 = 0.1, RCA_rca2 = 10, 
          Jmax25_rca2 = 50, Jmax25_npq4 = 50, TPU25_spsa = 3, Jmax25_rwt43 = 100)

parnames = names(pars)[1:13]


# Fit the model
# Don't use NPQ from rca2 for the fitting
aciPhoto = aci$Photo
rPhoto = transients$rPhoto
NPQ = transients$NPQ

Nobs = length(c(aciPhoto, rPhoto, NPQ[which(transients$NPQ > 0)]))

if(.Platform$OS.type == "windows") {
  cl <- makeCluster(6)
} else {
  cl <- makeForkCluster(6)
}
registerDoParallel(cl)
clusterEvalQ(cl, {library(MiniModel)})
clusterExport(cl, c("run_aci", "run_transient", "Ci", "inputs", "Nobs"))


fitAll <- nlfit(model = run_all, obs = c(aciPhoto, rPhoto, NPQ[which(transients$NPQ > 0)]),
                start = pars, algorithm = "subplex", fixed = "KiqEp",
                options = list(parscale = abs(pars), maxiter = 1e3, xrtol = 0, frtol = 1e-4),
                upper = upper, lower = lower, trace = "plot")

fitAll <- nlfit(model = run_all, obs = c(aciPhoto, rPhoto, NPQ[which(transients$NPQ > 0)]),
                   start = fitAll$all_par, algorithm = "subplex", fixed = "KiqEp",
                   options = list(parscale = abs(fitAll$all_par), maxiter = 1e3, xrtol = 0, frtol = 1e-4),
                   upper = upper, lower = lower, trace = "plot")



save(fitAll, file = "Intermediate/fitAll.RData")

load("Intermediate/fitAll.RData")


# Visualize fits
test = fitAll$model(fitAll$all_par)
aci = aci %>% mutate(simPhoto = test[1:nrow(aci)])
transients = transients %>% group_by() %>% mutate(simrPhoto = test[(nrow(aci) + 1):(nrow(aci) + nrow(transients))],
                                                  simNPQ = NPQ)
transients[which(transients$NPQ > 0),"simNPQ"][[1]] = test[(nrow(aci) + nrow(transients) + 1):length(test)]
transients[which(transients$TransientType == "transient1" & transients$Obs == 120),"NPQ"] = 1e-3
transients[which(transients$TransientType == "transient1" & transients$Obs == 120),"simNPQ"] = 1e-3

ggplot(aci, aes(x = Ci, y = Photo, colour = Genotype)) + geom_point() + geom_line(aes(y = simPhoto))
ggplot(transients, aes(x = Obs, y = rPhoto, colour = Genotype)) + geom_point() + facet_grid(.~TransientType) + geom_line(aes(y = simrPhoto))
ggplot(transients %>% filter(NPQ > 0), aes(x = Time - 3720, y = NPQ, colour = Genotype)) + 
  geom_point() + geom_line(aes(y = simNPQ)) + xlim(0,3600)


# Calculate residuals for each set
aci = aci %>% mutate(resPhoto = simPhoto - Photo)
transients = transients %>% mutate(resIS = simrPhoto - rPhoto, resNPQ = simNPQ - NPQ)

# Calculate statistics for each dataset
nAci = nrow(aci)
CDaci = with(aci, 1 - var(resPhoto)/var(Photo))
nRMSEaci = with(aci, var(resPhoto)/mean(Photo))

nIS = nrow(transients)
CDIS = with(transients, 1 - var(resIS)/var(rPhoto))
nRMSEIS = with(transients, var(resIS)/mean(rPhoto))

nNPQ = nrow(transients %>% filter(NPQ > 0))
CDNPQ = with(transients %>% filter(NPQ > 0), 1 - var(resNPQ)/var(NPQ))
nRMSENPQ = with(transients %>% filter(NPQ > 0), var(resNPQ)/mean(NPQ))

experiment_summary["ACI",] = c("all", nRMSEaci, CDaci, nAci)
experiment_summary["Induction",] = c("all", nRMSEIS, CDIS, nIS)
experiment_summary["NPQ",] = c("all", nRMSENPQ, CDNPQ, nNPQ)



# Helper function to create slice between two known values
# This is necessary as many parameters are one-side bounded
# The reason being that we impose a boundary on how much can NPQ increase at the end of the simulation
# To ensure that the model does not overestimate qI by "crossing" the data
slice2 = function(var, bounds, n = 20) {
  mid = fitAll$all_par[[var]]
  minNLL = fitAll$NLL
  steps = c(seq(bounds[1], mid, l = n/2 + 1)[-(n/2 + 1)],
            seq(mid, bounds[2], l = n/2 + 1)[-1])
  NLL = foreach(i = 1:n, .combine = "c") %do% {
    pars = fitAll$all_par
    pars[[var]] = steps[i]
    nlfit(model = run_all, obs = c(aciPhoto, rPhoto, NPQ[which(transients$NPQ > 0)]),
          start = pars, algorithm = "subplex", fixed = "KiqEp",
          options = list(parscale = abs(fitAll$all_par), maxiter = 1e3, xrtol = 0, frtol = 1e-4),
          upper = upper, lower = lower, eval = TRUE)$NLL
  }
  DNLL = NLL - minNLL
  profile = cbind(steps, DNLL)
  profile = cbind(profile, p = pchisq(q = profile[,2]*2, df = 1))
  profile = rbind(profile[1:(n/2),], cbind(mid, 0, 0), profile[(n/2 + 1):n,])
  profile
}

# Container for slices
sliceAll = vector("list")
class(sliceAll) = "profile"
attr(sliceAll, "level") = 0.95

# slice for KiqEp (no upper boundary)
sliceAll$KiqEp = slice2("KiqEp", c(0.05,0.12))

# slice for RB (no upper boundary)
sliceAll$RB = slice2("RB", fitAll$par["RB"]*c(0.995,1))

# slice for Jmax25
sliceAll$Jmax25 = slice2("Jmax25", fitAll$par["Jmax25"]*c(0.995,1.005))

# slice for TPU25
sliceAll$TPU25 = slice2("TPU25", fitAll$par["TPU25"]*c(0.995,1.005))

# slice for PhiqEmax (no lower boundary)
sliceAll$PhiqEmax = slice2("PhiqEmax", fitAll$par["PhiqEmax"]*c(1,1.005))

# slice for Kinh0 (no upper boundary)
sliceAll$Kinh0 = slice2("Kinh0", fitAll$par["Kinh0"]*c(0.995,1))

# # slice for fprot (no lower boundary)
sliceAll$fprot = slice2("fprot", fitAll$par["fprot"]*c(1,1.01))

# # slice for KiqEz (no upper boundary)
sliceAll$KiqEz = slice2("KiqEz", fitAll$par["KiqEz"]*c(0.98,1))

# # slice for Krca
sliceAll$Krca = slice2("Krca", fitAll$par["Krca"]*c(0.99,1.01))

# slice for fRBmin
sliceAll$fRBmin = slice2("fRBmin", fitAll$par["fRBmin"]*c(0.8,1.20))

# slice for Vrmax
sliceAll$Vrmax = slice2("Vrmax", fitAll$par["Vrmax"]*c(0.98,1.02))

# slice for gamma1
sliceAll$gamma1 = slice2("gamma1", fitAll$par["gamma1"]*c(0.98,1.02))

# slice for gamma2 (no upper boundary)
sliceAll$gamma2 = slice2("gamma2", fitAll$par["gamma2"]*c(0.98,1))

plot(sliceAll)


stopCluster(cl)

save(sliceAll, file = "Intermediate/sliceAll.RData")

# Some bounds cannot be calculated due to boundaries or flattened slices
# Others fail because the shape cannot be well represented by a spline
load("Intermediate/sliceAll.RData")
confintAll = confint(sliceAll)
confintAll["KiqEp", "upper"] = Inf
confintAll["KiqEp", "lower"] = approx(sliceAll$KiqEp[4:5,2], sliceAll$KiqEp[4:5,1], 2)$y
confintAll["RB", "upper"] = 2*confintAll["RB", "value"] - confintAll["RB", "lower"]
confintAll["PhiqEmax", "lower"] = 2*confintAll["PhiqEmax", "value"] - confintAll["PhiqEmax", "upper"]
confintAll["Kinh0", "upper"] = 2*confintAll["Kinh0", "value"] - confintAll["Kinh0", "lower"]
confintAll["fprot", "lower"] = 2*confintAll["fprot", "value"] - confintAll["fprot", "upper"]
confintAll["KiqEz", "upper"] = 2*confintAll["KiqEz", "value"] - confintAll["KiqEz", "lower"]
confintAll["gamma2", "upper"] = 2*confintAll["gamma2", "value"] - confintAll["gamma2", "lower"]

confintAll["RCA_rca2",] = c(fitAll$par["RCA_rca2"], NA, NA)
confintAll["Jmax25_rca2",] = c(fitAll$par["Jmax25_rca2"], NA, NA)
confintAll["Jmax25_npq4",] = c(fitAll$par["Jmax25_npq4"], NA, NA)
confintAll["Jmax25_rwt43",] = c(fitAll$par["Jmax25_rwt43"], NA, NA)
confintAll["TPU25_spsa",] = c(fitAll$par["TPU25_spsa"], NA, NA)

all(confintAll[,"lower"] < confintAll[,"value"], na.rm = T)
all(confintAll[,"upper"] > confintAll[,"value"], na.rm = T)
confintAll

confintAll["KdqEp",] = confintAll["KiqEp",]
confintAll["KdqEz",] = confintAll["KiqEz",]

for(i in row.names(confintAll))
  params_experiment[i,] = confintAll[i,]
params_experiment["RCA_rca2",] = c(unname(fitAll$par["RCA_rca2"]), NA, NA)
params_experiment["Jmax25_npq4",] = c(unname(fitAll$par["Jmax25_npq4"]), NA, NA)
params_experiment["Jmax25_rca2",] = c(unname(fitAll$par["Jmax25_rca2"]), NA, NA)
params_experiment["Jmax25_rwt43",] = c(unname(fitAll$par["Jmax25_rwt43"]), NA, NA)
params_experiment["TPU25_spsa",] = c(unname(fitAll$par["TPU25_spsa"]), NA, NA)

save(params_experiment, experiment_summary, file = "Intermediate/ParametersExperiment.RData")




