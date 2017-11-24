# This file assigns values to parameters, either by taking the values directly from the literature or fitting to some data

# Libraries ---------------------------------------------------------------
library(LeafLicor)
library(nlfit)
library(readr)
library(dplyr)
library(RcppSundials)

# Extract default parameters ----------------------------------------------
model = generate_LeafLicor_model()
default_params = as.list(model$Parameters$Values/model$Parameters$Coefs)
params_literature = data.frame(value = rep(NA, length(default_params)), lower = NA, upper = NA)
row.names(params_literature) = names(default_params)
data_summary = data.frame(params = NA, NRMSE = NA, CD = NA, n = NA)


# Value takes directly from literature or result of calculations ----------
params_literature["alphab", "value"] = 0.92
params_literature["alphabp", "value"] = 0.66
params_literature["alphag", "value"] = 0.72
params_literature["alphagp", "value"] = 0.60
params_literature["alphared", "value"] = 0.83
params_literature["alpharp", "value"] = 0.80
params_literature["falphaSc", "value"] = 0.93
params_literature["fcyc", "value"] = 0.05
params_literature["fpseudo", "value"] = 0.13
params_literature["gcm25", "value"] = 0.39
params_literature["gw25", "value"] = 0.75
params_literature["D0", "value"] = 1 # This is assumed, but reasonable (see Leuning, 1995)
params_literature["kD0", "value"] = 4.55e8
params_literature["kDinh", "value"] = 5e9
params_literature["kf", "value"] = 6.9e7
params_literature["KiPGA", "value"] = 0.84
params_literature["KmPGA", "value"] = 5
params_literature["KmRuBP", "value"] = 2.2e-2
params_literature["kp", "value"] = 2.65e9
params_literature["RCA", "value"] = 117.37
params_literature["Scm", "value"] = 7.1
params_literature["Sm", "value"] = 9.8
params_literature["sigma2", "value"] = 0.5
params_literature["leaf_surface", "value"] = 2
params_literature["theta", "value"] = 0.7
params_literature["Vch", "value"] = 7
params_literature["Vref", "value"] = 1.5e-4
params_literature["O2", "value"] = 210
params_literature["volume_chamber", "value"] = 80
params_literature["Flow", "value"] = 500
params_literature["gbw", "value"] = 9.2
params_literature["KdqEp", "value"] = 3.000000e-02
params_literature["KdqEz", "value"] = 4.830000e-04

# Temperature response of RCA ---------------------------------------------

# Data from Carmo-Silva and Salvucci (2011)
data = data_frame(Tleaf   = c(20, 22.43, 24.98, 27.49, 30, 32.54, 34.94, 37.53, 39.96) + 273.15,
                  fRCA = c(0.12, 0.14, 0.14, 0.15, 0.13, 0.11, 0.08, 0.04, 0.01))

# Fit peaked function to data
R = 8.31
start = c(RCA_m = 0.14, RCA_Hd = 1e5, Topt = 302, RCA_Ha = 1e4)
fit = nlfit(model = fRCA~RCA_m*RCA_Hd*exp((Tleaf - Topt)*RCA_Ha/(Topt*R*Tleaf))/(RCA_Hd - RCA_Ha*(1 - exp(RCA_Hd*(Tleaf - Topt)/(Topt*R*Tleaf)))),
            start = start, lower = c(0.10, 0,280,0), upper = c(0.20, Inf, 340, Inf),
            algorithm = "subplex", data = data, options = list(parscale = start))

# Calculate confidence intervals from profile likelihood
prof = profile(fit, parm = names(start), lower_delta = c(0.1, 0.1, 0.01, 0.2))
plot(prof)
confint_Temp_Rca = confint(prof)

# Assign parameters
params_literature["DHaRCA",] = confint_Temp_Rca["RCA_Ha",]/1e3
params_literature["DHdRCA",] = confint_Temp_Rca["RCA_Hd",]/1e3
params_literature["ToRCA",] = confint_Temp_Rca["Topt",]

# Write summary statistics from the fitting
data_summary["RCA_Temp",] = c(paste0(c("DHaRCA", "DHdRCA", "ToRCA"), collapse = ", "), 
                              NRMSE = round(rmse(fit)/mean(data$fRCA), 3), 
                              CD = round(cd(fit), 3), n = nrow(data))


# Temperature response of Kmc ---------------------------------------------

# Data from Walker et al (2013)
data = read_csv("Input/DigitizedData/Walker2013_Kmc.csv", col_names = c("Tleaf","Kmc"), col_types = "dd") %>%
          mutate(Kmc = Kmc/101e3*1e6, Tleaf = Tleaf + 273.15)

# Fit normalized Arrhenius
R = 8.31
Tref = 298.15
start = c(H = 35.4e3, Kmc25 = 200)
fit = nlfit(model = Kmc~Kmc25*exp((H*(Tleaf - Tref))/(Tref*R*Tleaf)),
            start = start, lower = c(0,0), algorithm = "subplex", 
            data = data, options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["Kmc25",] = conf["Kmc25",]
params_literature["DHaKmc",] = conf["H",]/1e3

# Write summary statistics from the fitting
data_summary["Kmc_Temp",] = c(paste0(c("Kmc25", "DHaKmc"), collapse = ", "), 
                              NRMSE = round(rmse(fit)/mean(data$Kmc), 3), 
                              CD = round(cd(fit), 3), n = nrow(data))


# Temperature response of Kmo ---------------------------------------------

# Data from Walker et al (2013)
data = read_csv("Input/DigitizedData/Walker2013_Kmo.csv", col_names = c("Tleaf","Kmo"), col_types = "dd") %>%
          mutate(Kmo = Kmo/101e3*1e6, Tleaf = Tleaf + 273.15)

# Fit normalized Arrhenius
R = 8.31
Tref = 298.15
start =  c(H = 35.4e3, Kmo25 = 200)
fit = nlfit(model = Kmo~Kmo25*exp((H*(Tleaf - Tref))/(Tref*R*Tleaf)),
            start = start, lower = c(0, 0), algorithm = "subplex", 
            data = data, options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["Kmo25",]  = conf["Kmo25",]
params_literature["DHaKmo",] = conf["H",]/1e3


# Write summary statistics from the fitting
data_summary["Kmo_Temp",] = c(paste0(c("Kmo25", "DHaKmo"), collapse = ", "), 
                              NRMSE = round(rmse(fit)/mean(data$Kmo), 3), 
                              CD = round(cd(fit), 3), n = nrow(data))



# Temperature response of Vcmax ---------------------------------------------

# Data from Walker et al (2013) on normalized Kc
data = read_csv("Input/DigitizedData/Walker2013_Vcmax.csv", col_names = c("Tleaf","Kc"), col_types = "dd") %>%
          mutate(Kc = Kc*4.4, Tleaf = Tleaf + 273.15)

# Fit normalized Arrhenius
R = 8.31
Tref = 298.15
start = c(H = 40e3, Kc25 = 1)
fit = nlfit(model = Kc~Kc25*exp((H*(Tleaf - Tref))/(Tref*R*Tleaf)),
            start = start, lower = c(0, 0), algorithm = "subplex", 
            data = data, options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["DHaKc",] = conf["H",]/1e3
params_literature["Kc25",] = conf["Kc25",]


# Write summary statistics from the fitting
data_summary["Kc_Temp",] = c(paste0(c("Kc25", "DHaKc"), collapse = ", "), 
                              NRMSE = round(rmse(fit)/mean(data$Kc), 3), 
                              CD = round(cd(fit), 3), n = nrow(data))


# Temperature response of Vomax ---------------------------------------------

# Data from Walker et al (2013) on normalized Ko
data = read_csv("Input/DigitizedData/Walker2013_Vomax.csv", col_names = c("Tleaf","Ko"), col_types = "dd") %>%
          mutate(Ko = Ko*4.4, Tleaf = Tleaf + 273.15)

# Fit normalized Arrhenius
R = 8.31
Tref = 298.15
start = c(H = 40e3, Ko25 = 0.8)
fit = nlfit(model = Ko~Ko25*exp((H*(Tleaf - Tref))/(Tref*R*Tleaf)),
            start = start, lower = c(0, 0), algorithm = "subplex", 
            data = data, options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["DHaKo",] = conf["H",]/1e3
params_literature["Ko25",] = params_literature["Kc25",]*0.2 # This is a more common Vomax/Vcmax ratio and it produces the right gamma_star


# Write summary statistics from the fitting
data_summary["Ko_Temp",] = c(paste0(c("Ko25", "DHaKo"), collapse = ", "), 
                             NRMSE = round(rmse(fit)/mean(data$Ko), 3), 
                             CD = round(cd(fit), 3), n = nrow(data))



# Temperature response of Rd ---------------------------------------------

# Data from Walker et al (2013) on normalized Ko
data = read_csv("Input/DigitizedData/Walker2013_Rd.csv", col_names = c("Tleaf","Rd"), col_types = "dd") %>% 
          mutate(Rd = Rd*4.4, Tleaf = Tleaf + 273.15)

# Fit normalized Arrhenius
R = 8.31
Tref = 298.15
start = c(H = 40e3, Rd25 = 1)
fit = nlfit(model = Rd~Rd25*exp((H*(Tleaf - Tref))/(Tref*R*Tleaf)),
            start = start, lower = c(0, 0), algorithm = "subplex", 
            data = data, options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["DHaRm",] = conf["H",]/1e3


# Write summary statistics from the fitting
data_summary["Rd_Temp",] = c(paste0(c("DHaRm"), collapse = ", "), 
                             NRMSE = round(rmse(fit)/mean(data$Rd), 3), 
                             CD = round(cd(fit), 3), n = nrow(data))


# Effect temperature on kinetics of chloroplast movement ------------------

# Data fromn Brugnoli et al. (1992)
Brugnoli = read_csv("Input/DigitizedData/Brugnoli.csv", col_names = c("Tleaf", "Rate"), col_types = "dd") %>% 
  mutate(Tleaf = Tleaf + 273.15)

# Fit peaked function to data
R = 8.31
Tref = 298.15
start = c(p25 = 4, Ha = 9e4, Topt = 305, Hd = 3.28e5)
fit = nlfit(model = Rate~p25*Hd*exp((Tleaf - Topt)*Ha/(Topt*R*Tleaf))/(Hd - Ha*(1 - exp(Hd*(Tleaf - Topt)/(Topt*R*Tleaf)))),
            start = start, lower = c(0, 0,300,0), upper = c(6, Inf, 308, Inf), 
            algorithm = "subplex", data = Brugnoli,
            options = list(maxiter = 1e6, xrtol = 0, frtol = 0, parscale = start))
Topt = fit$par[["Topt"]]
Ha = fit$par[["Ha"]]
Hd = fit$par[["Hd"]]
start = c(p25 = 4, Ha = Ha, Ds = Hd/Topt + R*log(Ha/(Hd - Ha)), Hd = Hd)
fit = nlfit(model = Rate~p25*exp((Tleaf - Tref)*Ha/(Tref*R*Tleaf))*
              (1 + exp((Tref*Ds - Hd)/(Tref*R)))/(1 + exp((Tleaf*Ds - Hd)/(Tleaf*R))),
            start =start, algorithm = "subplex", data = Brugnoli,
            options = list(parscale = start, maxiter = 1e6, xrtol = 0))
while(TRUE) {
  NLL = fit$NLL
  fit = update(fit)
  if(abs((NLL - fit$NLL)/NLL) < 1e-10) break
}

prof = profile(fit, parm = names(start), lower_delta = c(0.01,0.1, 0.05, 0.1), algorithm = "praxis",
               options = list(randseed = 1L, xrtol = 1e-10))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["DHaKalpha",] = conf["Ha",]/1e3
params_literature["DHdKalpha",] = conf["Hd",]/1e3
params_literature["DsKalpha",] = conf["Ds",]/1e3


# Write summary statistics from the fitting
data_summary["Chl_Temp",] = c(paste0(c("DHaKalpha", "DHdKalpha", "DsKalpha"), collapse = ", "), 
                             NRMSE = round(rmse(fit)/mean(Brugnoli$Rate), 3), 
                             CD = round(cd(fit), 3), n = nrow(Brugnoli))


# Steady-state and kinetics of chloroplast movement -----------------------

# Data from Davis et al. (2012) and Labuz et al. (2015)
Davis2012 = read_csv(file = "Input/DigitizedData/Davis2012.csv", col_names = c("Blue","DeltaTau"), col_types = "dd") %>%
  mutate(alphar = 1 - DeltaTau*1.3826/100)
Labuz2015_23 = read_csv("Input/DigitizedData/Labuz2015_23degrees.csv", col_names = c("Time","DeltaTau"), col_types = "dd") %>%
  arrange(Time) %>% mutate(Stage = c(rep(0,3),rep(1,18),rep(2,14),rep(3,18)),
                           alphar = 1 - DeltaTau*1.3826/100,Time = Time*60,
                           Blue = c(rep(0,3),rep(1.6,18),rep(20,14),rep(120,18)))
Labuz2015_30 = read_csv("Input/DigitizedData/Labuz2015_30degrees.csv", col_names = c("Time","DeltaTau"), col_types = "dd") %>%
  arrange(Time) %>% mutate(Stage = c(rep(0,3),rep(1,13),rep(2,14),rep(3,17)),
                           alphar = 1 - DeltaTau*1.3826/100,Time = Time*60,
                           Blue = c(rep(0,3),rep(1.6,13),rep(20,14),rep(120,17)))  
Labuz2015_4 = read_csv("Input/DigitizedData/Labuz2015_4degrees.csv", col_names = c("Time","DeltaTau"), col_types = "dd") %>%
  arrange(Time) %>% mutate(Stage = c(rep(0,4),rep(1,18),rep(2,24),rep(3,13)),
                           alphar = 1 - DeltaTau*1.3826/100,Time = Time*60,
                           Blue = c(rep(0,4),rep(1.6,18),rep(20,24),rep(120,13))) 

# Model of steady-state relative absorptance of red light as a function of blue light and temperature
steady_alphar <- function(parms, Ib, T) {
  with(as.list(parms), {
    R = 8.31
    Tref = 298.15
    alphar_alpha = alphar_alpha25*exp(-(T - Tref)*DHaAlphar/(Tref*R*T))
    alpharss = pmin(1.0 + Ib/Iac*alpharac, 
                   1.0 + alpharac - (alphar_alpha*(Ib - Iac) + alpharav - sqrt((alphar_alpha*(Ib - Iac) + alpharav)^2 -
                      4*alphar_alpha*thetaalphar*alpharav*(Ib - Iac)))/(2*thetaalphar))
    if(any(is.nan(alpharss))) stop("Negative sqrt")
    alpharss
  })
}

setwd("Code/Calibration/ChloroplastModel")
try(dyn.unload("ChloroplastModel"))
system("R CMD SHLIB ChloroplastModel.cpp")
dyn.load("ChloroplastModel")
setwd("../../..")
dynamic_alphar = getNativeSymbolInfo(name = "chloroplast")$address

# Objective function to fit steady-state and dynamic data
model = function(parms) {
  
  # NLL for steady-state data
  steady = try(steady_alphar(parms, Davis2012$Blue, 298.15), silent = TRUE)
  if(inherits(steady, "try-error")) return(1e12)
  
  # NLL for dynamic data (use a different amplitude for light avoidance movement)
  parms["alpharav"] = parms["alpharav"]*parms["falpharav"]
  
  settings =  SimulationModels:::create_settings_ODE(list(atol = 1e-12, rtol = 1e-12), 1, 0, 13)
  
  # Simulate alphar for each temperature, using lsoda from deSolve and returning a very large number if there is an error in the simulation
  sim23 = try(wrap_cvodes(times = Labuz2015_23$Time, states_ = c(alphar = 1), parameters_ =  c(parms[-13], T = 273.15 + 23), 
                      forcings_data_ = list(Blue23), model_ = dynamic_alphar, settings = settings, jacobian_ = dynamic_alphar))
  if(inherits(sim23, "try-error")) return(1e12)

  sim4 = try(wrap_cvodes(times = Labuz2015_4$Time, states_ = c(alphar = 1), parameters_ =  c(parms[-13], T = 273.15 + 4), 
                     forcings_data_ = list(Blue4), model_ = dynamic_alphar, settings = settings, jacobian_ = dynamic_alphar))
  if(inherits(sim4, "try-error")) return(1e12)
  
  sim30 = try(wrap_cvodes(times = Labuz2015_30$Time, states_ = c(alphar = 1), parameters_ =  c(parms[-13], T = 273.15 + 30), 
                     forcings_data_ = list(Blue30), model_ = dynamic_alphar, settings = settings, jacobian_ = dynamic_alphar))
  if(inherits(sim30, "try-error")) return(1e12)
  
  dynamic = c(sim23[,2], sim4[,2], sim30[,2])
  
  c(dynamic, steady)
  
}

# Default parameter values (we have run the optimization several times)
parms = c(alphar_alpha25 = 6.556156e-03, 
          DHaAlphar = 9.151325e+04, 
          Iac = 2.347008e+00, 
          alpharac = 5.346126e-02 , 
          alpharav = 2.501486e-01, 
          thetaalphar = 1e-6, 
          alphar_K_i = 2.244587e-03, 
          alphar_K_d = 2.895236e-03, 
          HaKalpha = params_literature["DHaKalpha","value"]*1e3, 
          DsKalpha = params_literature["DsKalpha","value"]*1e3, 
          HdKalpha = params_literature["DHdKalpha","value"]*1e3, 
          falpharav = 1.065050)
lower = rep(0, length(parms))
upper = rep(Inf, length(parms))
names(lower) = names(parms)
names(upper) = names(parms)
upper[c(1, 4:6)] = 1
upper["Iac"] = 5
upper["falpharav"] = 2

ineq = function(parms) {
  if(any(is.nan(parms))) {return(1e12)}
  steady = try(steady_alphar(parms, Davis2012$Blue, 298.15), silent = TRUE)
  if(inherits(steady, "try-error")) {return(1e12)} else return(-1)
}

# Interpolators for light in each experiment
switches = Labuz2015_23[which(diff(Labuz2015_23$Blue) > 0) + 1,"Time"][[1]]
#Blue23 = approxfun(x = sort(c(0,switches, switches - 1,5e4)),y =c(0,0,1.6,1.6,20,20,120,120),method = "linear", rule = 2)
Blue23 = cbind(x = sort(c(0,switches, switches - 1, max(Labuz2015_23$Time))),y =c(0,0,1.6,1.6,20,20,120, 120))
switches = Labuz2015_4[which(diff(Labuz2015_4$Blue) > 0) + 1,"Time"][[1]]
#Blue4 = approxfun(x = sort(c(0,switches, switches - 1,5e4)),y =c(0,0,1.6,1.6,20,20,120,120),method = "linear", rule = 2)
Blue4 = cbind(x = sort(c(0,switches, switches - 1, max(Labuz2015_4$Time))),y =c(0,0,1.6,1.6,20,20,120,120))
switches = Labuz2015_30[which(diff(Labuz2015_30$Blue) > 0) + 1,"Time"][[1]]
#Blue30 = approxfun(x = sort(c(0,switches, switches - 1,5e4)),y =c(0,0,1.6,1.6,20,20,120,120),method = "linear", rule = 2)
Blue30 = cbind(x = sort(c(0,switches, switches - 1, max(Labuz2015_30$Time))),y =c(0,0,1.6,1.6,20,20,120,120))

# Test the model in C++
settings =  SimulationModels:::create_settings_ODE(list(), 1, 0, length(parms) + 1)
test = wrap_cvodes(times = 1:5000, states_ = c(alphar = 1), parameters_ =  c(parms, T = 273.15 + 23), forcings_data_ = list(as.matrix(Labuz2015_23[,c("Time", "Blue")])),
                   model_ = dynamic_alphar, settings = settings, jacobian_ = dynamic_alphar)

#Fit the model
# fit = nlfit(model = model, start = parms, fixed = c("HaKalpha", "DsKalpha", "HdKalpha"),
#                    lower = lower, upper = upper, ineq = ineq,
#                    algorithm = "subplex", obs = c(Labuz2015_23$alphar, Labuz2015_4$alphar, Labuz2015_30$alphar, Davis2012$alphar),
#                    options = list(parscale = parms, maxiter = 1e4, frtol = 1e-6, xrtol = 0))
# fit$lower["thetaalphar"] = -Inf
# # Difficulties to converge. Need to reset the algorithm several times
# for(i in 1:100) {
#   print(i)
#   NLL = fit$NLL
#   fit = update(fit, options = list(frtol = 1e-10, xrtol = 0))
#   if(abs(fit$NLL - NLL)/abs(NLL) < 1e-10) break
# }
# save(fit, file = "Intermediate/FitChloroplast.RData")

load("Intermediate/FitChloroplast.RData")
# Use Hessian approximation (n = 170)
conf = confint(fit, parm = c("alphar_alpha25", "DHaAlphar", "Iac", "alpharac", "alpharav", 
                               "thetaalphar", "alphar_K_i", "alphar_K_d", "falpharav"), type = "hessian")


# Assign parameters to the model
params_literature["alphar_alpha25",] = conf["alphar_alpha25",]
params_literature["DHaAlphar",] = conf["DHaAlphar",]/1e3
params_literature["Iac",] = conf["Iac",]
params_literature["alpharac",] = conf["alpharac",]
params_literature["alpharav",] = conf["alpharav",]
params_literature["thetaalphar",] = conf["thetaalphar",]
params_literature["Kialpha25",] = conf["alphar_K_i",]
params_literature["Kdalpha25",] = conf["alphar_K_d",]


# Write summary statistics from the fitting
meas = c(Labuz2015_23$alphar, Labuz2015_4$alphar, Labuz2015_30$alphar, Davis2012$alphar)
data_summary["Chl_Dynamic",] = c(paste0(c("alphar_alpha25", "DHaAlphar", "Iac", "Kialpha25",
                                       "alpharac", "alpharav", "thetaalphar", "Kdalpha25"), collapse = ", "), 
                              NRMSE = round(rmse(fit)/mean(meas), 3), 
                              CD = round(cd(fit), 3), n = length(meas))

n1 = nrow(Labuz2015_23)
n2 = nrow(Labuz2015_4)
n3 = nrow(Labuz2015_30)
n4 = nrow(Davis2012)
sim23 = cbind(Labuz2015_23$Time, fit$fitted[1:n1])
sim4 = cbind(Labuz2015_4$Time, fit$fitted[(n1 + 1):(n1 + n2)])
sim30 = cbind(Labuz2015_30$Time, fit$fitted[(n1 + n2 + 1):(n1 + n2 + n3)])
simSteady = cbind(Davis2012$Blue, fit$fitted[(n1 + n2 + n3 + 1):(n1 + n2 + n3 + n4)])

windowsFonts(Cambria=windowsFont("Cambria Math"))
png("Output/ChloroplastMovement.png",, width = 24, height = 12, units = "cm", 
    res = 600, bg = "white", antialias = "default", family ="Cambria")

par(yaxs = "i", xaxs = "i", las = 1, cex.axis = 1.5, cex.lab = 1.5, lwd = 2, mgp = c(3.2,1,0),
    mar = c(4.8,4.8,1,0), mfrow = c(1,2))

plot(sim23[,1]/60,sim23[,2], ylim = c(0.80,1.10), yaxt = "n", xaxt = "n",
     xlim = c(-20,420), t = "l", ylab = expression(italic(alpha[r])),
     xlab = "Time (min)", lwd = 2)
points(Labuz2015_23$Time/60, Labuz2015_23$alphar)
lines(sim4[,1]/60,sim4[,2], col = 2, lty = 2, lwd = 2)
points(Labuz2015_4$Time/60, Labuz2015_4$alphar, col = 2, pch = 2)
lines(sim30[,1]/60,sim30[,2], col = 3, lty = 3, lwd = 2)
points(Labuz2015_30$Time/60, Labuz2015_30$alphar, col = 3, pch = 3)

axis(side = 1, at = seq(0,400,100), tcl = -0.5, lwd.ticks = 2)
axis(side = 1, at = seq(50,400,100), labels = NA, tcl = -0.25, lwd.ticks = 2)
axis(side = 2, at = seq(0.80,1.1,0.05), tcl = -0.5, lwd.ticks = 2)
axis(side = 2, at = seq(0.825,1.075,0.05), labels = NA, tcl = -0.25, lwd.ticks = 2)
axis(side = 3, at = seq(0,400,100), labels = NA, tcl = 0.5, lwd.ticks = 2)
axis(side = 3, at = seq(50,400,100), labels = NA, tcl = 0.25, lwd.ticks = 2)

legend("bottomleft", legend = c(expression(4~degree*C),
        expression(23~degree*C), expression(30~degree*C)),
       col = c(2,1,3), lty = c(2,1,3), pch = c(2,1,3), bty = "n", cex = 1.5)
text(x = 0, y = 1.075, labels = "A", cex = 1.5)

par(mar = c(4.8,0,1,4.8))
plot(Davis2012[,c("Blue","alphar")], ylim = c(0.80,1.1), xlim = c(-5,150),
     ylab = "", yaxt = "n",xaxt = "n",
     xlab = expression(italic(I[b])~(mu*mol~m^{-2}~s^{-1})))
lines(0:150, steady_alphar(fit$all_par, 0:150, 25 + 273.15), lwd = 2)

axis(side = 1, at = seq(0,150,30), lwd.ticks = 2, tcl = -.5)
axis(side = 1, at = seq(15,150,30), labels = NA, lwd.ticks = 2, tcl = -0.25)
axis(side = 3, at = seq(0,150,30), labels = NA, lwd.ticks = 2, tcl = .5)
axis(side = 3, at = seq(15,150,30), labels = NA, lwd.ticks = 2, tcl = 0.25)
axis(side = 4, at = seq(0.80,1.1,0.05), labels = NA, tcl = 0.5, lwd.ticks = 2)
axis(side = 4, at = seq(0.825,1.075,0.05), labels = NA, tcl = 0.25, lwd.ticks = 2)

text(x = 0, y = 1.075, labels = "B", cex = 1.5)

dev.off()



# Effect of light intensity on enzyme activity ----------------------------

# Data from Sassenrath-Cole et al (1994) on FBPase activity as a function of light intensity
SassenrathFBPase = read_csv("Input/DigitizedData/Sassenrath_FBPase.csv", col_names = c("I0","Vr"), col_types = "dd")

# Fit non-rectangular hyperbolae to the data
start = c(fR_alpha = 0.0025, fR_theta = 0.99, fR0 = 0.1, Vrmax = 24)
fit = nlfit(model = Vr~Vrmax*(fR0 + (fR_alpha*I0  + (1 - fR0) - 
                    sqrt((fR_alpha*I0 + (1 - fR0))^2 - 4*fR_alpha*fR_theta*I0*(1 - fR0)))/(2*fR_theta)), 
            start = start, algorithm = "subplex", 
            data = SassenrathFBPase,lower = rep(0, length(start)), upper = c(1, 1, 1, 100), 
            options = list(xrtol = 1e-10, parscale = start))
plot(fit, type = "response", var = SassenrathFBPase$I0)
fit$lower[3]= -1
prof = profile(fit, parm = names(start)[-2])
prof_theta = profile(fit, parm ="fR_theta", lower_delta = 0.01, direction = "lower")
prof$fR_theta = prof_theta[[1]]
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["fR0",] = conf["fR0",]
params_literature["alphafR",] = conf["fR_alpha",]
params_literature["thetafR",] = conf["fR_theta",]

data_summary["FBPase_Light",] = c(paste0(c("fR0", "alphafR", "thetafR"), collapse = ", "), 
                                 NRMSE = round(rmse(fit)/mean(SassenrathFBPase$Vr), 3), 
                                 CD = round(cd(fit), 3), n = nrow(SassenrathFBPase))

# Dynamics of enzyme activation -------------------------------------------

# Data from Sassenrath-Cole et al (1994) on dynamic changes in enzyme activity
SassenrathFBPaseUp = read_csv("Input/DigitizedData/Sassenrath_FBPase_time_up.csv", col_names = c("Time", "Vr"), col_types = "dd") %>%
                       arrange(Time) %>% mutate(fR = Vr/mean(Vr[7:8]))
SassenrathFBPaseDown = read_csv("Input/DigitizedData/Sassenrath_FBPase_time_down.csv", col_names = c("Time", "Vr"), col_types = "dd") %>%
                        arrange(Time) %>% mutate(fR = Vr/Vr[1]) 
SassenrathRubiscoDown = read_csv("Input/DigitizedData/Sassenrath_Rubisco_time_down.csv", col_names = c("Time", "Vc"), col_types = "dd") %>%
  arrange(Time) %>% mutate(fRB = Vc/Vc[1]) 

# Fit exponential transients to the data
start = c(k = 1/5/60, FBPasemin = 0.2)
fit = nlfit(fR~1 - (1 - FBPasemin)*exp(-Time*60*k), data = SassenrathFBPaseUp, 
            start = start, algorithm = "subplex",
            lower = c(0,0), upper = c(300, 1), options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)
params_literature["KiR",] = conf["k",]


data_summary["FBPase_Up",] = c(paste0(c("KiR"), collapse = ", "), 
                                  NRMSE = round(rmse(fit)/mean(SassenrathFBPaseUp$fR), 3), 
                                  CD = round(cd(fit), 3), n = nrow(SassenrathFBPaseUp))

fit = nlfit(fR~FBPasemin + (1 - FBPasemin)*exp(-Time*60*k), data = SassenrathFBPaseDown, 
                     start = start, algorithm = "subplex",
                     lower = c(0,0), upper = c(300, 1), options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)
params_literature["KdR",] = conf["k",]

data_summary["FBPase_Down",] = c(paste0(c("KdR"), collapse = ", "), 
                               NRMSE = round(rmse(fit)/mean(SassenrathFBPaseDown$fR), 3), 
                               CD = round(cd(fit), 3), n = nrow(SassenrathFBPaseDown))


start = c(k = 1/3600, RBmin = 0.2)
fit = nlfit(fRB~RBmin + (1 - RBmin)*exp(-Time*60*k), data = SassenrathRubiscoDown, 
            start = start, algorithm = "subplex",
            lower = c(k = 0, RBmin = 0.1), upper = c(k = 1, RBmin = 0.3), options = list(parscale = start))
fit$lower["RBmin"] = -1
plot(fit, type = "response", var = SassenrathRubiscoDown$Time*60)
prof = profile(fit, parm = "k")
plot(prof)
conf = confint(prof)
params_literature["KdRB",] = conf["k",]

data_summary["Rubisco_Down",] = c(paste0(c("KdRB"), collapse = ", "), 
                               NRMSE = round(rmse(fit)/mean(SassenrathRubiscoDown$fRB), 3), 
                               CD = round(cd(fit), 3), n = nrow(SassenrathRubiscoDown))


# Rca effect on steady-state Rubisco activation ---------------------------

# Data from Mott and Woodrow (2000)
Mott2000RCA = read_csv("Input/DigitizedData/Mott_2000_KaRCA.csv", col_names = c("RCA","An"), col_types = "dd")

# Fit Michaelis-Menten equation to relationship between steady-state Rubisco activity and Rca (using RB-limited An as surrogate)
start = c(KaRCA = 14, Anmax = 20)
fit = nlfit(An~Anmax*RCA/(RCA + KaRCA), data = Mott2000RCA, 
            start = start, algorithm = "subplex",
            lower = c(0,0), upper = c(30, 100), options = list(parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters to model
params_literature["KaRCA",] = conf["KaRCA",]

data_summary["KaRCA",] = c(paste0(c("KaRCA"), collapse = ", "), 
                                  NRMSE = round(rmse(fit)/mean(Mott2000RCA$An), 3), 
                                  CD = round(cd(fit), 3), n = nrow(Mott2000RCA))

# Rate constant of PSII repair --------------------------------------------

Kasahara = read_csv("Input/DigitizedData/Kasahara.csv", col_names = c("Time", "FvFm"), col_types = "dd") %>%
  arrange(Time) %>% mutate(Time = (Time - Time[1])*3600)

# At 70 umol/m2/s, phiqE is almost 0 so there is no photoprotection, so kinh = kinh0 (around 0.16)
parms = c(fIIa0 = 0.3, kr = 1e-4)

model = function(pars) {
  fIIa0 = pars[[1]]
  krep = pars[[2]]
  fIIa = 1 - (1 - fIIa0)*exp(-krep*Kasahara$Time)
  kf = params_literature["kf","value"]
  kD0 = params_literature["kD0","value"]
  kDinh = params_literature["kDinh","value"]
  kp = params_literature["kp","value"]
  Fm  = kf/(kf + kD0)
  Fmd = kf/(kf + kDinh)
  Fo  = kf/(kf + kD0 + kp)
  FvFm = (Fm - Fo)*fIIa/(Fm*fIIa + Fmd*(1 - fIIa))
}
# Fit analytical model of photoinhibition
fit = nlfit(model = model, obs = Kasahara$FvFm, start = parms, algorithm = "subplex",
            options = list(parscale = parms))
prof = profile(fit, parm = names(parms))
plot(prof)
conf = confint(prof)


# Assign parameters to model
params_literature["Krep25",] = conf["krep",]

data_summary["Krep",] = c(paste0(c("krep"), collapse = ", "), 
                           NRMSE = round(rmse(fit)/mean(Kasahara$FvFm), 3), 
                           CD = round(cd(fit), 3), n = nrow(Kasahara))


# Rate constant photorespiration ------------------------------------------

# Time series of glycine decrease in the darkness normalized by initial value
rawsthorne = data.frame(Time = c(0,30.45,44.40, 278.26),
                          Glycine = c(0.95,0.60,0.33, 0.13))
rawsthorne$Glycine = rawsthorne$Glycine/0.95

# Fit rate constant
start = c(kpr = 1/50, G0 = 0.1)
fit = nlfit(model = Glycine~G0 + (1 - G0)*exp(-Time*kpr), data = rawsthorne, start = start, 
            algorithm = "subplex", options = list(xrtol = 1e-10, parscale = start))
prof = profile(fit, parm = names(start))
plot(prof)
conf = confint(prof)

# Assign parameters to model
params_literature["kPR",] = conf["kpr",]

data_summary["kPR",] = c(paste0(c("kPR"), collapse = ", "), 
                          NRMSE = round(rmse(fit)/mean(rawsthorne$Glycine), 3), 
                          CD = round(cd(fit), 3), n = nrow(rawsthorne))


# Temperature response of Jmax --------------------------------------------

# Data from Walker et al (2013)
data = read_csv("Input/DigitizedData/Yamori_JmaxHT.csv", col_names = c("Temperature", "Jmax"), col_types = "dd") %>%
                mutate(Temperature = Temperature + 273.15)


# Fit normalized modified Arrhenius
R = 8.31
Tref = 298.15
start = c(Jmax25 = 110, Jmax_Ha = 1e4, Jmax_Ds = 1e2, Jmax_Hd = 1e5)
fit = nlfit(model = Jmax~Jmax25*exp((Temperature - Tref)*Jmax_Ha/(Tref*R*Temperature))*
              (1 + exp((Tref*Jmax_Ds - Jmax_Hd)/(Tref*R)))/(1 + exp((Temperature*Jmax_Ds - Jmax_Hd)/(Temperature*R))),
            start = start, algorithm = "subplex", data = data,
            lower = rep(0, length(start)),
            options = list(parscale = start, xrtol = 0, xatol = 0, fatol = 0, frtol = 1e-10, 
                           maxiter = 1e6))
prof = profile(fit, parm = names(start), lower_delta = c(0.2, 0.2, 0.6, 0.2))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["Jmax25",]  = conf["Jmax25",]
params_literature["DHaJmax",] = conf["Jmax_Ha",]/1e3
params_literature["DsJmax",] = conf["Jmax_Ds",]/1e3
params_literature["DHdJmax",] = conf["Jmax_Hd",]/1e3

data_summary["Jmax_Temp",] = c(paste0(c("Jmax25", "DHaJmax", "DsJmax", "DHdJmax"), collapse = ", "), 
                         NRMSE = round(rmse(fit)/mean(data$Jmax), 3), 
                         CD = round(cd(fit), 3), n = nrow(data))


# Temperature response of TPU ---------------------------------------------

# Data from Harley et al (1992)
data = read_csv("Input/DigitizedData/HarleyTPU.csv", col_names = c("Tleaf", "TPU"), col_types = "dd") %>%
  mutate(Tleaf = Tleaf + 273.15) %>% arrange(Tleaf)


# Fit Arrhenius with optimal temperature (easier than original one) and the calculate the original model
R = 8.31
Tref = 298.15
start = c(TPU_o = 11, TPU_Ha = 6e4, Topt = 302, TPU_Hd = 2e5)
fit = nlfit(model = TPU~TPU_o*TPU_Hd*exp((Tleaf - Topt)*TPU_Ha/(Topt*R*Tleaf))/(TPU_Hd - TPU_Ha*(1 - exp(TPU_Hd*(Tleaf - Topt)/(Topt*R*Tleaf)))),
            start = start, algorithm = "subplex", data = data,
            lower = c(10,1e4, 300, 1e5), upper = c(12,1e5,306,5e5),
            options = list(parscale = start, maxiter = 1e6, xrtol = 1e-10))
Topt = fit$par[["Topt"]]
TPU_Ha = fit$par[["TPU_Ha"]]
TPU_Hd = fit$par[["TPU_Hd"]]
start = c(TPU25 = 7.5, TPU_Ha = TPU_Ha, TPU_Ds = TPU_Hd/Topt + R*log(TPU_Ha/(TPU_Hd - TPU_Ha)), TPU_Hd = TPU_Hd)
fit = nlfit(model = TPU~TPU25*exp((Tleaf - Tref)*TPU_Ha/(Tref*R*Tleaf))*
              (1 + exp((Tref*TPU_Ds - TPU_Hd)/(Tref*R)))/(1 + exp((Tleaf*TPU_Ds - TPU_Hd)/(Tleaf*R))),
            start =start, algorithm = "subplex", data = data,
            lower = c(7,0.95*start[-1]), upper = c(9,1.05*start[-1]),
            options = list(parscale = start, maxiter = 1e6, xrtol = 1e-10))
plot(data)
lines(290:310, predict(fit, newdata = data.frame(Tleaf = 290:310)))
fit$lower = rep(0, length(fit$all_par)); names(fit$lower) = names(fit$all_par);
fit$upper = rep(Inf, length(fit$all_par)); names(fit$upper) = names(fit$all_par)
prof = nlfit:::slicer(fit, parm = names(start), lower_delta = c(0.2, 0.2, 0.01, 0.01), nstep = 30)

conf = confint(prof)


# Assign parameters
params_literature["TPU25",]  = conf["TPU25",]
params_literature["DHaTPU",] = conf["TPU_Ha",]/1e3
params_literature["DsTPU",] = conf["TPU_Ds",]/1e3
params_literature["DHdTPU",] = conf["TPU_Hd",]/1e3


data_summary["TPU_Temp",] = c(paste0(c("TPU25", "DHaTPU", "DsTPU", "DHdTPU"), collapse = ", "), 
                               NRMSE = round(rmse(fit)/mean(data$TPU), 3), 
                               CD = round(cd(fit), 3), n = nrow(data))



# Temperature response of PSII repair -------------------------------------

# Data from Greer et al (1986)
data = read_csv("Input/DigitizedData/GreerKrep.csv", col_names = c("Tleaf", "Krep"), col_types = "dd") %>%
  mutate(Tleaf = Tleaf + 273.15, Krep = Krep*1e-2/60)

# Fit normalized modified Arrhenius
R = 8.31
Tref = 298.15
start = c(Krep25 = 1e-4, Krep_Ha = 60e3, Krep_Ds = 0.6e3, Krep_Hd = 201e3)
fit = nlfit(model = Krep~Krep25*exp((Tleaf - Tref)*Krep_Ha/(Tref*R*Tleaf))*(1 + exp((Tref*Krep_Ds - Krep_Hd)/(Tref*R)))/
              (1 + exp((Tleaf*Krep_Ds - Krep_Hd)/(Tleaf*R))),
            start = start, algorithm = "subplex", data = data,
            options = list(parscale = start, xrtol = 0, maxiter = 1e6))
plot(fit, type = "response", var = data$Tleaf)
prof = profile(fit, parm = names(start), options = list(maxiter = 1e4, xrtol = 1e-10, dx = 0.1), lower_delta = c(0.2, 0.2, 0.05,0.2))
plot(prof)
conf = confint(prof)

# Assign parameters
params_literature["Krep25",]  = conf["Krep25",]
params_literature["DHaKrep",] = conf["Krep_Ha",]/1e3
params_literature["DsKrep",] = conf["Krep_Ds",]/1e3
params_literature["DHdKrep",] = conf["Krep_Hd",]/1e3


data_summary["Krep_Temp",] = c(paste0(c("Krep25", "DHaKrep", "DsKrep", "DHdKrep"), collapse = ", "), 
                              NRMSE = round(rmse(fit)/mean(data$Krep), 3), 
                              CD = round(cd(fit), 3), n = nrow(data))


# Temperature response gm -------------------------------------------------

# Data from von Caemmerer (2015)
data = read_csv("Input/DigitizedData/gmTemp.csv", col_names = c("Tleaf", "gm"), col_types = "dd") %>%
  mutate(Tleaf = Tleaf + 273.15, gm = gm*0.13)

# Fit normalized modified Arrhenius
R = 8.31
Tref = 298.15
start = c(gm25 = 0.13, gm_Ha = 66.3e3, gm_Ds = 0.3e3, gm_Hd = 91e3)
fit = nlfit(model = gm~gm25*exp((Tleaf - Tref)*gm_Ha/(Tref*R*Tleaf))*(1 + exp((Tref*gm_Ds - gm_Hd)/(Tref*R)))/
              (1 + exp((Tleaf*gm_Ds - gm_Hd)/(Tleaf*R))),
            start = start, algorithm = "subplex", data = data,
            options = list(parscale = start, xrtol = 0, frtol = 0, maxiter = 1e6))
plot(fit, type = "response", var = data$Tleaf)
prof = profile(fit, parm = names(start), options = list(maxiter = 1e4, xrtol = 1e-10, dx = 0.05, randseed = 1), 
               lower_delta = c(0.2, 0.2, 0.3,0.3), algorithm = "praxis")
plot(prof)
conf = confint(prof)

par(las = 1)
plot(data, ylab = "gm", xlab = "Tleaf", ylim = c(0,0.15),xlim = c(285,310))
Tl = seq(285,310,0.1)
lines(Tl, with(as.list(fit$par), gm25*exp((Tl - Tref)*gm_Ha/(Tref*R*Tl))*(1 + exp((Tref*gm_Ds - gm_Hd)/(Tref*R)))/
  (1 + exp((Tl*gm_Ds - gm_Hd)/(Tl*R)))))


# Assign parameters
# params_literature["gcm25",]  = conf["gm25",]
params_literature["DHaGc",] = conf["gm_Ha",]/1e3
params_literature["DsGc",] = conf["gm_Ds",]/1e3
params_literature["DHdGc",] = conf["gm_Hd",]/1e3
# params_literature["gw25",]  = conf["gm25",]
params_literature["DHaGw",] = conf["gm_Ha",]/1e3
params_literature["DsGw",] = conf["gm_Ds",]/1e3
params_literature["DHdGw",] = conf["gm_Hd",]/1e3


data_summary["gm_Temp",] = c(paste0(c("gcm25", "DHaGc", "DsGc", "DHdGc",
                                        "gcw25", "DHaGw", "DsGw", "DHdGw"), collapse = ", "), 
                               NRMSE = round(rmse(fit)/mean(data$gm), 3), 
                               CD = round(cd(fit), 3), n = nrow(data))

# Effect of VPD on gs -----------------------------------------------------
# In the new version, the value is assumed (the data was not related to At anyways)

# # Data from Dai et al (1992)
# data = read_csv("Input/DigitizedData/DaiGs.csv", col_names = c("VPD","gs"), col_types = "dd")
# 
# # Fit hyperbola
# start = c(D0 = 0.5, gsmax = 400)
# fit = nlfit(model = gs~gsmax/(1 + VPD/D0),
#             start = start, algorithm = "subplex", data = data,
#             options = list(parscale = start, xrtol = 0, maxiter = 1e6))
# plot(fit, type = "response", var = data$VPD)
# prof = profile(fit, parm = names(start))
# plot(prof)
# conf = confint(prof)
# 
# # Assign parameters
# params_literature["D0",]  = conf["D0",]
# 
# 
# data_summary["gs_VPD",] = c(paste0(c("D0"), collapse = ", "), 
#                              NRMSE = round(rmse(fit)/mean(data$gs), 3), 
#                              CD = round(cd(fit), 3), n = nrow(data))


# Fit effect of CO2 on Rubisco activity -----------------------------------

# Data from von Caemmerer et al. (1986)
# Cc/Ci ratio from Caemmerer and Evans (1991) on the same species
data = read_csv("Input/DigitizedData/Caemmerer1986_1.csv", col_names = c("Ci","fECM"), col_types = "dd") %>%
          mutate(fECM = fECM/100, Cc = Ci*0.65) 

# Fit linear model to low CO2 data
fit = lm(fECM~Cc, data = subset(data, Cc < 60))
conf = confint(fit)

params_literature["ac",] = c(coef(fit)[1], conf[1,])
params_literature["bc",] = c(coef(fit)[2], conf[2,])


data_summary["Carbamylation",] = c(paste0(c("ac", "bc"), collapse = ", "), 
                            NRMSE = round(rmse(fit)/mean(data$fECM), 3), 
                            CD = round(summary(fit)$adj.r.squared, 3), n = nrow(data))



# Save data to external file ----------------------------------------------
save(params_literature, data_summary, file = "Intermediate/ParametersLiterature.RData")

