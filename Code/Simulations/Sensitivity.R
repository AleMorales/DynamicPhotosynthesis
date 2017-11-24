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


# Fluctuating light -----------------------------------------------------------------------------------------------


lightflecks = function(model, period, PARs) {
  
  PAR1 = PARs[1]
  PAR2 = PARs[2]
  
  # CVODE settings
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                     c(1e-14,1e-6,1e4,20,20,20, -1e-6))
  model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  
  # Assign parameters from literature to the model
  filter = which(names(parameters) %in% names(model$Parameters$Values))
  model$set_parameters(names(parameters[filter]),unname(parameters[filter])) 
  
  # Simulate a LiCOR with 10% red - 90% blue
  model$set_forcings("Ib", cbind(c(0,1), c(5,5)))
  model$set_forcings("Ir", cbind(c(0,1), c(45,45)))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("CO2R", cbind(c(0,1), c(400,400)))
  model$set_forcings("H2OR", cbind(c(0,1), c(20,20)))
  model$set_forcings("Ta", cbind(c(0,1), c(298.15, 298.15)))
  model$set_forcings("Tl", cbind(c(0,1), c(298.15, 298.15)))
  model$set_states("Ci", 400)
  model$set_states("Cc", 400)
  model$set_states("Ccyt", 400)
  model$set_states("PR", 25)
  
  # Calculate steady-state
  model$set_time(c(0,3600))
  steadyState = cvode(model)[2,names(model$States$Values)]
  model$set_states(names(steadyState), steadyState)
  
  # Turn on sensitivity analysis
  model$set_settings("sensitivity", TRUE)
  model$set_settings("which_sens", which(names(model$Parameters$Values) %in% 
                       c("KiR", "KdR", "Krca", "KdRB", "Kgsi", "Kgsd", "Kinh0",
                         "Krep25", "KiqEp", "KdqEp", "KiqEz", "KdqEz", "Kialpha25",
                         "Kdalpha25")))
  model$set_settings("which_observed", which(model$Observed$Names == "A"))
  model$set_settings("which_states", which(names(model$States$Values) %in% c("gsw","Ca","Ci")))
  model$set_settings("maxtime", 600)
  
  # Simulate Square wave light determined by period
  timeI = sort(c(seq(0,3600, l = 3600/period), seq(1e-2,3600 + 1e-2, l = 3600/period)))
  model$set_forcings("Ib", cbind(timeI, rep(c(PAR1,PAR2,PAR2,PAR1)*0.1, times = 3600/period/2)))
  model$set_forcings("Ir", cbind(timeI, rep(c(PAR1,PAR2,PAR2,PAR1)*0.9, times = 3600/period/2)))
  model$set_time(seq(0,3600,by = min(period/10, 1)))
  sim = cvode(model)
  class(sim) = "matrix"
  sim = as.data.frame(sim)
  
  # Calculate sensitive of A to each rate constant
  S_A_KiqEp = with(sim, (S_Ca_KiqEp - S_Ci_KiqEp)*gsw + (Ca - Ci)*S_gsw_KiqEp)
  S_A_KdqEp = with(sim, (S_Ca_KdqEp - S_Ci_KdqEp)*gsw + (Ca - Ci)*S_gsw_KdqEp)
  S_A_KiqEz = with(sim, (S_Ca_KiqEz - S_Ci_KiqEz)*gsw + (Ca - Ci)*S_gsw_KiqEz)
  S_A_KdqEz = with(sim, (S_Ca_KdqEz - S_Ci_KdqEz)*gsw + (Ca - Ci)*S_gsw_KdqEz)
  S_A_Kinh0 = with(sim, (S_Ca_Kinh0 - S_Ci_Kinh0)*gsw + (Ca - Ci)*S_gsw_Kinh0)
  S_A_Krep25 = with(sim, (S_Ca_Krep25 - S_Ci_Krep25)*gsw + (Ca - Ci)*S_gsw_Krep25)
  S_A_Kialpha25 = with(sim, (S_Ca_Kialpha25 - S_Ci_Kialpha25)*gsw + (Ca - Ci)*S_gsw_Kialpha25)
  S_A_Kdalpha25 = with(sim, (S_Ca_Kdalpha25 - S_Ci_Kdalpha25)*gsw + (Ca - Ci)*S_gsw_Kdalpha25)
  S_A_KiR = with(sim, (S_Ca_KiR - S_Ci_KiR)*gsw + (Ca - Ci)*S_gsw_KiR)
  S_A_KdR = with(sim, (S_Ca_KdR - S_Ci_KdR)*gsw + (Ca - Ci)*S_gsw_KdR)
  S_A_Krca = with(sim, (S_Ca_Krca - S_Ci_Krca)*gsw + (Ca - Ci)*S_gsw_Krca)
  S_A_KdRB = with(sim, (S_Ca_KdRB - S_Ci_KdRB)*gsw + (Ca - Ci)*S_gsw_KdRB)
  S_A_Kgsi = with(sim, (S_Ca_Kgsi - S_Ci_Kgsi)*gsw + (Ca - Ci)*S_gsw_Kgsi)
  S_A_Kgsd = with(sim, (S_Ca_Kgsd - S_Ci_Kgsd)*gsw + (Ca - Ci)*S_gsw_Kgsd)
  return(data.frame(S_A_KiqEp = S_A_KiqEp, S_A_KdqEp = S_A_KdqEp, S_A_KiqEz = S_A_KiqEz, 
                    S_A_KdqEz = S_A_KdqEz, S_A_Kinh0 = S_A_Kinh0, S_A_Krep25 = S_A_Krep25,
                    S_A_Kialpha25 = S_A_Kialpha25, S_A_Kdalpha25 = S_A_Kdalpha25,
                    S_A_KiR = S_A_KiR, S_A_KdR = S_A_KdR, S_A_Krca = S_A_Krca, S_A_KdRB = S_A_KdRB,
                    S_A_Kgsi = S_A_Kgsi, S_A_Kgsd = S_A_Kgsd, A = sim[,'A']))
}

# Simulate all mutants and periods -------------------------------------------------------------------------------------

# if(.Platform$OS.type == "windows") {
#   cl <- makeCluster(8) # Hyperthreading does not offset overhead of sockets
# } else {
#   cl <- makeForkCluster(8)
# }
# registerDoParallel(cl)
# 
periods = c(0.1, 0.5, 1, 3, 5, 15, 30, 60, 120, 360)
# LFsens = foreach(period = periods, .packages = "MiniModel") %dopar% {
#   lightflecks(generate_MiniModel_model(), period, c(50, 1000))
# }
# 
# stopCluster(cl)
# 
# save(LFsens, file = "Intermediate/LFsens.RData")
# 


load("Intermediate/LFsens.RData")

# Compute relative change in A per relative in input

senspars = parameters[c("KiqEp", "KdqEp",  "KiqEz",  "KdqEz",  "Kinh0", "Krep25", 
                        "Kialpha25",  "Kdalpha25", "KiR",  "KdR",  "Krca",  "KdRB", 
                        "Kgsi", "Kgsd")]
W = length(senspars)

# Calculate increase in A per relative change in input (units of dA)
LFDA = plyr::llply(LFsens, function(x) {
  L = nrow(x)
  out = cbind(x[,W + 1], x[,1:W]*matrix(rep(senspars, L), nrow = L, ncol = W, byrow = T))
  colnames(out) = c("A",names(senspars))
  out
})

# Then the integrate A and dA over time
LFDAtotal = plyr::laply(LFDA, colSums)

# Now calculate the relative change in A per relative change in input
LFeffect = LFDAtotal/LFDAtotal[,1]

LFeffect = cbind(Period = periods, LFeffect)

png("Output/figureLF.png", width = 7, height = 12, pointsize = 8, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(2,1), xaxs = "i", yaxs = "i", las = 1, mar= c(0.2,0,0,0), oma = c(4,4,0.5,1))

with(as.data.frame(LFeffect), {
  plot(Period, KiR, t = "l", log = "x", ylim = c(-0.02,0.55), xaxt = "n", xlim = c(0.1,365),
       xlab = "Period (s)", ylab  = expression(Delta*A))
  lines(Period, Krca, col = 2, lty = 2, t = "l")
  lines(Period, Kgsi, col = 3, lty = 3, t = "l")
})
legend("topright", c("KiR", "RB", "Kgsi"), col = 1:3, lty = 1:3, bty = "n", ncol = 1)

with(as.data.frame(LFeffect), {
  plot(Period, KdqEp, t = "l", log = "x", ylim = c(-0.001,0.015), xaxt = "n", xlim = c(0.1,365),
       xlab = "", ylab  = expression(Delta*A), yaxt = "n")
  lines(Period, KdqEz, col = 2, lty = 2, t = "l")
  lines(Period, Krep25, col = 3, lty = 3, t = "l")
  lines(Period, Kialpha25, col = 4, lty = 4, t = "l")
  axis(1, at = c(0.1,1,5,30,120,360), labels = c("0.1","1","5","30","120","360"), cex.axis = 0.9)
  axis(2, seq(0, 0.015, 0.005))
  legend("top", c("KdqEp", "KdqEz",  "Krep25", "Kdalpha25"), col = 1:4, lty = 1:4, bty = "n", ncol = 2)
})
mtext("Period (s)", side = 1, line = 2.5, outer = T)

dev.off()
