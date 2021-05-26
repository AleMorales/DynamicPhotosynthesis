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


# Induction & Relaxation curve ---------------------------------------------------------

run = function(model) {
  
  # CVODE settings
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                     c(1e-10,1e-6,1e4,20,20,20, -1e-6))
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
  tryCatch(model$set_states("PR", 25), error = function(x) NULL)
  
  # Calculate steady-state
  model$set_time(c(0,3600))
  steadyState = cvode(model)[2,names(model$States$Values)]
  model$set_states(names(steadyState), steadyState)
  
  # Simulate Transient
  model$set_forcings("Ib", cbind(c(0,3600, 3601, 7200), c(100,100,5,5)))
  model$set_forcings("Ir", cbind(c(0,3600, 3601, 7200), c(900,900,45,45)))
  model$set_time(1:7200)
  Transient = cvode(model)
  class(Transient) = "matrix"
  as_data_frame(Transient)
}



# Simulate all mutants --------------------------------------------------------------------------------------------

control = run(generate_MiniModel_model())
QssfR = run(generate_MiniModelQssfR_model())
QssfRB = run(generate_MiniModelQssfRB_model())
QssKgs = run(generate_MiniModelQssgs_model())
QssqE = run(generate_MiniModelQssqE_model())
QssqM = run(generate_MiniModelQssqM_model())
QssqI = run(generate_MiniModelQssqI_model())
QssPR = run(generate_MiniModelQssPR_model())
Qss = run(generate_MiniModelQss_model())

# Calculate differences in Photo ----------------------------------------------------------------------------------

mutations = data_frame(fR = (QssfR$A - control$A),#/control$A,
                       fRB = (QssfRB$A - control$A),#/control$A,
                       Kgs = (QssKgs$A - control$A),#/control$A,
                       qE = (QssqE$A - control$A),#/control$A,
                       qI = (QssqI$A - control$A),#/control$A,
                       qM = (QssqM$A - control$A),#/control$A)
                       PR = (QssPR$A - control$A),#/control$A)
                       QSS = (Qss$A - control$A))#/control$A)



filter1 = 2:3600
filter2 = 3615:7200
time1 = (1:7200)[filter1]/60
time2 = (1:7200)[filter2]/60

png("Output/figureLimitations.png", width = 10, height = 6, pointsize = 8, units = "cm", 
    res = 600, bg = "white", antialias = "default")

with(mutations, {
  par(mfrow = c(1,1), xaxs = "i", yaxs = "i", las = 1, mar= c(4.0,4.0,0.5,1), mgp = c(2,1,0))
  
  plot(1,1, xlim = c(0,120), ylim = c(-0.3,5), type = "n", ylab = expression(italic(Delta*A)~(mu*mol~m^{-2}~s^{-1})),
       xlab = "Time (min)")
  lines(time1, fR[filter1],  col = 1, lty = 1)
  lines(time1, fRB[filter1], col = 2, lty = 2)
  lines(time1, Kgs[filter1], col = 3, lty = 3)
  lines(time2, fR[filter2],  col = 1, lty = 1)
  lines(time2, fRB[filter2], col = 2, lty = 2)
  lines(time2, Kgs[filter2], col = 3, lty = 3)
  lines(time1, qE[filter1],  col = 4, lty = 4)
  lines(time1, qI[filter1],  col = 5, lty = 5)
  lines(time1, qM[filter1],  col = 6, lty = 6)
  lines(time2, qE[filter2],  col = 4, lty = 4)
  lines(time2, qI[filter2],  col = 5, lty = 5)
  lines(time2, qM[filter2],  col = 6, lty = 6)
  lines(time1, PR[filter1],  col = 7, lty = 7)
  lines(time2, PR[filter2],  col = 7, lty = 7)
  lines(time1, QSS[filter1], col = 8, lty = 8)
  lines(time2, QSS[filter2], col = 8, lty = 8)
  abline(v = 60, lty = 2)
  abline(h = 0, lty = 2)
  legend("topright", c("R", "RB", "gs", "qE", "qI", "qM", "PR" ,"QSS"), col = 1:8, lty = 1:8, bty = "n")
  text(30, 4.8, "1000", cex = 1.3)
  text(90, 4.8, "50", cex = 1.3)
})

dev.off()


png("Output/figureCumulativeLimitations.png", width = 10, height = 6, pointsize = 8, units = "cm", 
    res = 600, bg = "white", antialias = "default")

with(mutations, {
  par(mfrow = c(1,1), xaxs = "i", yaxs = "i", las = 1, mar= c(4.0,5.2,0.5,1), mgp = c(2.7,1,0))
  
  plot(1,1, xlim = c(0,120), ylim = c(-100,3000), type = "n", ylab = expression(sum(italic(Delta*A)["i"]~(mu*mol~m^{-2}), i == 0, i == t)),
       xlab = "Time (min)")
  lines(time1, cumsum(fR[filter1]),  col = 1, lty = 1)
  lines(time1, cumsum(fRB[filter1]), col = 2, lty = 2)
  lines(time1, cumsum(Kgs[filter1]), col = 3, lty = 3)
  lines(time2, cumsum(fR[filter2]),  col = 1, lty = 1)
  lines(time2, cumsum(fRB[filter2]), col = 2, lty = 2)
  lines(time2, cumsum(Kgs[filter2]), col = 3, lty = 3)
  lines(time1, cumsum(qE[filter1]),  col = 4, lty = 4)
  lines(time1, cumsum(qI[filter1]),  col = 5, lty = 5)
  lines(time1, cumsum(qM[filter1]),  col = 6, lty = 6)
  lines(time2, cumsum(qE[filter2]),  col = 4, lty = 4)
  lines(time2, cumsum(qI[filter2]),  col = 5, lty = 5)
  lines(time2, cumsum(qM[filter2]),  col = 6, lty = 6)
  lines(time1, cumsum(PR[filter1]),  col = 7, lty = 7)
  lines(time2, cumsum(PR[filter2]),  col = 7, lty = 7)
  lines(time1, cumsum(QSS[filter1]), col = 8, lty = 8)
  lines(time2, cumsum(QSS[filter2]), col = 8, lty = 8)
  abline(v = 60, lty = 2)
  abline(h = 0, lty = 2)
  legend("topright", c("R", "RB", "gs", "qE", "qI", "qM", "PR", "QSS"), col = 1:8, lty = 1:8, bty = "n")
  text(30, 2800, "1000", cex = 1.3)
  text(90, 2800, "50", cex = 1.3)
})

dev.off()


# Fluctuating light -----------------------------------------------------------------------------------------------


lightflecks = function(model, period, PARs, param = NULL, variable = "A") {
  
  PAR1 = PARs[1]
  PAR2 = PARs[2]
  
  # CVODE settings
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                     c(1e-14,1e-8,1e4,20,20,20, -1e-6))
  model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  model$set_settings("maxtime", 1000)
  
  # Assign parameters from literature to the model
  filter = which(names(parameters) %in% names(model$Parameters$Values))
  model$set_parameters(names(parameters[filter]),unname(parameters[filter])) 
  
  # Assign parameter for sensitivity analysis
  if(!is.null(param)) model$set_parameters(names(param), param)
  
  # Simulate a LiCOR with 10% red - 90% blue
  model$set_forcings("Ib", cbind(c(0,1), c(PAR1,PAR1)*0.1))
  model$set_forcings("Ir", cbind(c(0,1), c(PAR1,PAR1)*0.9))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("CO2R", cbind(c(0,1), c(400,400)))
  model$set_forcings("H2OR", cbind(c(0,1), c(20,20)))
  model$set_forcings("Ta", cbind(c(0,1), c(298.15, 298.15)))
  model$set_forcings("Tl", cbind(c(0,1), c(298.15, 298.15)))
  model$set_states("Ci", 400)
  model$set_states("Cc", 400)
  model$set_states("Ccyt", 400)
  tryCatch(model$set_states("PR", 25), error = function(x) NULL)
  
  # Calculate steady-state
  model$set_time(c(0,1800))
  steadyState = cvode(model)[2,names(model$States$Values)]
  model$set_states(names(steadyState), steadyState)
  
  # Simulate Transient - Square wave light determined by period
  timeI = sort(c(seq(0,1800, by = period), seq(1e-2,1800 + 1e-2, by = period)))
  model$set_forcings("Ib", cbind(timeI, c(rep(c(PAR1,PAR2,PAR2,PAR1)*0.1, times = 1800/period/2), PAR1*0.1,PAR2*0.1)))
  model$set_forcings("Ir", cbind(timeI,  c(rep(c(PAR1,PAR2,PAR2,PAR1)*0.9, times = 1800/period/2), PAR1*0.9,PAR2*0.9)))
  model$set_time(seq(0,1800,by = min(period/10, 1)))
  cvode(model)[,variable]
}

# Calculate steady-state at 1000 and 50 uE

run_steady = function(model, PAR) {

  # CVODE settings
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"),
                     c(1e-14,1e-6,1e4,20,20,20, -1e-6))
  model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  model$set_settings("maxtime", 600)

  # Assign parameters from literature to the model
  filter = which(names(parameters) %in% names(model$Parameters$Values))
  model$set_parameters(names(parameters[filter]),unname(parameters[filter]))

  # Simulate a LiCOR with 10% red - 90% blue
  model$set_forcings("Ib", cbind(c(0,1), c(PAR,PAR)*0.1))
  model$set_forcings("Ir", cbind(c(0,1), c(PAR,PAR)*0.9))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("CO2R", cbind(c(0,1), c(400,400)))
  model$set_forcings("H2OR", cbind(c(0,1), c(20,20)))
  model$set_forcings("Ta", cbind(c(0,1), c(298.15, 298.15)))
  model$set_forcings("Tl", cbind(c(0,1), c(298.15, 298.15)))
  model$set_states("Ci", 400)
  model$set_states("Cc", 400)
  model$set_states("Ccyt", 400)
  tryCatch(model$set_states("PR", 25), error = function(x) NULL)

  # Calculate steady-state
  model$set_time(c(0,7200))
  steadyState = cvode(model)[2,]
}

test = run_steady(generate_MiniModel_model(), 0)
APARdata = matrix(NA, nrow = 30, ncol = length(test))
colnames(APARdata) = names(test)
PARdata = seq(0,1200, l = 30)
for(i in 1:30) {
  APARdata[i,] = run_steady(generate_MiniModel_model(), PARdata[i])
}
APAR = loess(A~PAR, data = data.frame(A = APARdata[,"A"], PAR = PARdata), span = 0.5,
             control = loess.control(surface = "direct", statistics = "exact"))
plot(PARdata, APARdata[,"A"])
lines(PARdata, predict(APAR, newdata=data.frame(PAR = PARdata)))

# Calculate QSS curve with all the rate constants increased
test = run_steady(generate_MiniModelQss_model(), 0)
APARdata2 = matrix(NA, nrow = 30, ncol = length(test))
colnames(APARdata2) = names(test)
model = generate_MiniModel_model()
filter = which(names(parameters) %in% names(model$Parameters$Values))
model$set_parameters(names(parameters[filter]),unname(parameters[filter]))
for(i in 1:30) {
  APARdata2[i,] = run_steady(generate_MiniModelQss_model(), PARdata[i])
}
plot(PARdata, APARdata2[,"A"])
points(PARdata, APARdata[,"A"], col = 2)
APAR = loess(A~PAR, data = data.frame(A = APARdata2[,"A"], PAR = PARdata), span = 0.5,
             control = loess.control(surface = "direct", statistics = "exact"))
RpPAR = loess(A~PAR, data = data.frame(A = APARdata2[,"Rp"], PAR = PARdata), span = 0.5,
             control = loess.control(surface = "direct", statistics = "exact"))
RuBPPAR = loess(A~PAR, data = data.frame(A = APARdata2[,"RuBP"], PAR = PARdata), span = 0.5,
             control = loess.control(surface = "direct", statistics = "exact"))
fRBPAR = loess(A~PAR, data = data.frame(A = APARdata2[,"fRB"], PAR = PARdata), span = 0.5,
                control = loess.control(surface = "direct", statistics = "exact"))

lines(PARdata, predict(APAR, newdata=data.frame(PAR = PARdata)))

# Simulate all mutants and periods -------------------------------------------------------------------------------------

if(.Platform$OS.type == "windows") {
  cl <- makeCluster(8)
} else {
  cl <- makeForkCluster(8)
}
registerDoParallel(cl)


pars = list(parameters["KiR"],
            parameters["Krca"],
            parameters["Kgsi"],
            parameters[c("KdqEp", "KdqEz")],
            parameters["Krep25"],
            parameters["Kialpha25"],
            parameters["kPR"])
W = length(pars)


periods = c(0.1, 0.3, 0.5, 1, 2, 3,  5, 10, 15, 30, 60, 90, 150, 225)

LF = foreach(period = periods, .packages = "MiniModel") %dopar% {
  N = length(seq(0,1800,by = min(period/10, 1)))
  out = matrix(NA, ncol = W + 4, nrow = N)
  colnames(out) = c("control", "R", "RB", "gs", "qE", "qI", "qM", "PR", "All", "QSS", "QSS2")
  PARdata = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "PAR")
  out[,1] = lightflecks(generate_MiniModel_model(), period, c(50,1000))
  for(i in 2:(W + 1)) out[,i] = lightflecks(generate_MiniModel_model(), period, c(50,1000), pars[[i - 1]]*1e6)
  out[,W + 2] = lightflecks(generate_MiniModel_model(), period, c(50,1000), unlist(pars)*1e6)
  out[,W + 3] = lightflecks(generate_MiniModelQss_model(), period, c(50,1000))
  out[,W + 4] = PARdata
  out
}



stopCluster(cl)

# Check that average PAR is always the same
for(i in 1:length(LF)) {
  print(mean(LF[[i]][,"QSS2"]))
}


for(i in 1:length(LF)) {
  LF[[i]][,"QSS2"] = predict(APAR, newdata=data.frame(PAR = as.numeric(LF[[i]][,"QSS2"])))
}

save(LF, file = "Intermediate/LFTotal.RData")

load("Intermediate/LFTotal.RData")



# Compute relative change in A due to virtual mutant
processLF = function(LF) {
  LFeffect = LF
  LFhigh = LF
  LFlow = LF
  for(i in 1:length(LFeffect)) {
    period = periods[i]
    time = sort(c(seq(0,1800, l = 1800/period), seq(1e-2,1800 + 1e-2, l = 1800/period)))
    PAR = rep(c(0,1,1,0), times = 1800/period/2)
    PAR = approx(time, PAR, seq(0,1800,by = min(period/10, 1)))$y
    print(c(length(PAR), nrow(LFeffect[[i]])))
    # Time-integrated CO2 assimilation
    LFeffect[[i]] = colSums(LFeffect[[i]])
    LFeffect[[i]] = (LFeffect[[i]] - LFeffect[[i]][1])/LFeffect[[i]][1]
    
    # Low irradiance average CO2 assimilation
    LFlow[[i]] = colSums(LFlow[[i]][which(PAR <= 0.5),])
    LFlow[[i]] = (LFlow[[i]] - LFlow[[i]][1])/LFlow[[i]][1]
    
    # High irradiance average CO2 assimilation
    LFhigh[[i]] = colSums(LFhigh[[i]][which(PAR >= 0.5),])
    LFhigh[[i]] = (LFhigh[[i]] - LFhigh[[i]][1])/LFhigh[[i]][1]
  }
  LFeffect = do.call("rbind", LFeffect)
  LFeffect = cbind(Period = periods, LFeffect)
  colnames(LFeffect)[1] = "Period"
  
  LFlow = do.call("rbind", LFlow)
  LFlow = cbind(Period = periods, LFlow)
  colnames(LFlow)[1] = "Period"
  
  LFhigh = do.call("rbind", LFhigh)
  LFhigh = cbind(Period = periods, LFhigh)
  colnames(LFhigh)[1] = "Period"
  
  return(list(LFeffect, LFhigh, LFlow))
}

LFlist = processLF(LF)
LFeffect = LFlist[[1]]
LFhigh = LFlist[[2]]
LFlow = LFlist[[3]]

png("Output/figureLF.png", width = 14, height = 7, pointsize = 10, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(1,2), xaxs = "i", yaxs = "i", las = 1, mar = c(4.0,4.2,0.5,1), mgp = c(2,1,0))

with(as.data.frame(LFeffect), {
  plot(Period, R*100, t = "l", log = "x", ylim = c(-2,40), xaxt = "n", xlim = c(0.1,300),
       col = 1, lty = 1, 
       xlab = "Lightfleck duration (s)", ylab  = expression(italic(Delta*A/A)~("%")))
  lines(Period, RB*100, col = 2, lty = 2, t = "l")
  lines(Period, gs*100, col = 3, lty = 3, t = "l")
  lines(Period, qE*100, col = 4, lty = 4, t = "l")
  lines(Period, qI*100, col = 5, lty = 5, t = "l")
  lines(Period, qM*100, col = 6, lty = 6, t = "l")
  lines(Period, PR*100, col = 7, lty = 7, t = "l")
  axis(1, at = c(0.1, 1, 10, 100,200,300), labels = c("0.1","1", "10","100","", "300"))
  axis(1, at = 300, labels = "300")
  axis(1, at = seq(0.2,0.9,0.1), tcl = -0.2, labels = NA)
  axis(1, at = seq(2,9,1), tcl = -0.2, labels = NA)
  axis(1, at = seq(20,100,10), tcl = -0.2, labels = NA)
  legend("topleft", c("R", "RB", "gs", "qE", "qI", "qM", "PR"), col = 1:7, lty = 1:7, 
         bty = "n", ncol = 2, cex = 0.65, x.intersp = 0.5)
  abline(h= 0, lty = 1, col = "gray")
  text(10,38,labels = "A", cex = 1.2)
})

with(as.data.frame(LFeffect), {
  plot(Period, All*100, t = "l", log = "x", ylim = c(-10,40), xaxt = "n", xlim = c(0.1,300),
       col = 1, lty = 1, 
       xlab = "Lightfleck duration (s)", ylab  = expression(italic(Delta*A/A)~("%")))
  lines(Period, QSS2*100, col = 2 , lty = 2, t = "l")
  #lines(Period, QSS2*100, col = 3 , lty = 3, t = "l")
  axis(1, at = c(0.1, 1, 10, 100,200,300), labels = c("0.1","1", "10","100","", "300"))
  axis(1, at = 300, labels = "300")
  axis(1, at = seq(0.2,0.9,0.1), tcl = -0.2, labels = NA)
  axis(1, at = seq(2,9,1), tcl = -0.2, labels = NA)
  axis(1, at = seq(20,100,10), tcl = -0.2, labels = NA)
  legend("topleft", c("All", "QSS"), col = 1:2, lty = 1:2, 
         bty = "n", ncol = 1, cex = 0.65, x.intersp = 0.5)
  abline(h= 0, lty = 1, col = "gray")
  text(10,38,labels = "B", cex = 1.2)
})


dev.off()






# Plot dynamics of RuBP and photosynthesis during lightflecks -------------

# Get dynamics of RuBP
pars = list(parameters["KiR"],
            parameters["Krca"],
            parameters["Kgsi"],
            parameters[c("KdqEp", "KdqEz")],
            parameters["Krep25"],
            parameters["Kialpha25"],
            parameters["kPR"])
W = length(pars)
period = 0.5
N = length(seq(0,1800,by = min(period/10, 1)))

# RuBP dynamics
RuBP = matrix(NA, ncol = 3, nrow = N)
colnames(RuBP) = c("All", "control", "QSS2")
PARdata = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "PAR")
RuBP[,1] = lightflecks(generate_MiniModel_model(), period, c(50,1000), unlist(pars)*1e6, variable = "RuBP")
RuBP[,2] = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "RuBP")
RuBP[,3] = PARdata
RuBP[,3] = predict(RuBPPAR, newdata=data.frame(PAR = as.numeric(RuBP[,3])))

period = 90
N = length(seq(0,1800,by = min(period/10, 1)))
RuBP2 = matrix(NA, ncol = 3, nrow = N)
colnames(RuBP2) = c("All", "control", "QSS2")
PARdata = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "PAR")
RuBP2[,1] = lightflecks(generate_MiniModel_model(), period, c(50,1000), unlist(pars)*1e6, variable = "RuBP")
RuBP2[,2] = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "RuBP")
RuBP2[,3] = PARdata
RuBP2[,3] = predict(RuBPPAR, newdata=data.frame(PAR = as.numeric(RuBP2[,3])))

# Rubisco activity
period = 0.5
N = length(seq(0,1800,by = min(period/10, 1)))
fRB = matrix(NA, ncol = 3, nrow = N)
colnames(fRB) = c("All", "control", "QSS2")
PARdata = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "PAR")
fRB[,1] = lightflecks(generate_MiniModel_model(), period, c(50,1000), unlist(pars)*1e6, variable = "fRB")
fRB[,2] = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "fRB")
fRB[,3] = PARdata
fRB[,3] = predict(fRBPAR, newdata=data.frame(PAR = as.numeric(fRB[,3])))
PARdata5 = PARdata

period = 90
N = length(seq(0,1800,by = min(period/10, 1)))
fRB2 = matrix(NA, ncol = 3, nrow = N)
PARdata = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "PAR")
colnames(fRB2) = c("All", "control", "QSS2")
fRB2[,1] = lightflecks(generate_MiniModel_model(), period, c(50,1000), unlist(pars)*1e6, variable = "fRB")
fRB2[,2] = lightflecks(generate_MiniModel_model(), period, c(50,1000), variable = "fRB")
fRB2[,3] = PARdata
fRB2[,3] = predict(fRBPAR, newdata=data.frame(PAR = as.numeric(fRB2[,3])))




# Supplemental figure S6

png("Output/figureLightFlecksSample.png", width = 14, height = 18, pointsize = 10, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfcol = c(3,2), mar = c(4,4.5,0.5,0.65), las = 1, xaxs = "i", yaxs = "i")

t0 = 10000 + 11
t1 = t0 + 20
l = t0:t1 - t0

# Period= 0.5 s
plot(0.05*l, LF[[3]][t0:t1,"All"], t = "l", ylim = c(0,14),
     ylab = expression(A~(mu*mol~m^{-2}~s^{-1})), 
     xlab = "Time (s)")
lines(0.05*l, LF[[3]][t0:t1,"control"], col = 2)
minQSS2 = min(LF[[3]][t0:t1,"QSS2"])
maxQSS2 = max(LF[[3]][t0:t1,"QSS2"])
QSStime = c(0,0,10,10,20)
lines(0.05*QSStime, c(maxQSS2, minQSS2, minQSS2, maxQSS2, maxQSS2), col = 3)
abline(v = 0.05*(which(abs(diff(PARdata5[t0:t1])) > 100) - 1), lty = 2)
text(0.05,13.5, "A")
text(c(0.25, 0.75),13.5, c(50,1000))

legend("bottomright", c("control", "All", "QSS"), col = c(2,1,3), lty = 1, bty = "n")

plot(0.05*l, RuBP[t0:t1,"All"], t = "l", ylim = c(0,100),
     ylab = expression(RuBP~(mu*mol~m^{-2}~s^{-1})), 
     xlab = "Time (s)")
lines(0.05*l, RuBP[t0:t1,"control"], col = 2)
minQSS2 = min(RuBP[t0:t1,"QSS2"])
maxQSS2 = max(RuBP[t0:t1,"QSS2"])
QSStime = c(0,0,10,10,20)
lines(0.05*QSStime, c(maxQSS2, minQSS2, minQSS2, maxQSS2, maxQSS2), col = 3)
abline(v = 0.05*(which(abs(diff(PARdata5[t0:t1])) > 100) - 1), lty = 2)
text(0.05,95, "B")
text(c(0.25, 0.75),95, c(50,1000))


plot(0.05*l, fRB[t0:t1,"All"], t = "l", ylim = c(0,1),
     ylab = expression(fRB~(mu*mol~m^{-2}~s^{-1})), 
     xlab = "Time (s)")
lines(0.05*l, fRB[t0:t1,"control"], col = 2)
minQSS2 = min(fRB[t0:t1,"QSS2"])
maxQSS2 = max(fRB[t0:t1,"QSS2"])
QSStime = c(0,0,10,10,20)
lines(0.05*QSStime, c(maxQSS2, minQSS2, minQSS2, maxQSS2, maxQSS2), col = 3)
abline(v = 0.05*(which(abs(diff(PARdata5[t0:t1])) > 100) - 1), lty = 2)
text(0.05,0.95, "C")
text(c(0.25, 0.75),0.95, c(50,1000))


#  Period = 90 s
t0 = 1801 - 270
t1 = 1801 - 90
l = t0:t1 - t0

plot(l, LF[[12]][t0:t1,"All"], t = "l", ylim = c(0,14),
     ylab = expression(A~(mu*mol~m^{-2}~s^{-1})), 
     xlab = "Time (s)", xaxt = "n")
axis(1, seq(0,180,30))
lines(l, LF[[12]][t0:t1,"control"], col = 2)
minQSS2 = min(LF[[12]][t0:t1,"QSS2"])
maxQSS2 = max(LF[[12]][t0:t1,"QSS2"])
QSStime = c(0,0,90,90,180)
lines(QSStime, c(maxQSS2, minQSS2, minQSS2, maxQSS2, maxQSS2), col = 3)
abline(v = (which(abs(diff(PARdata[t0:t1])) > 100) - 1), lty = 2)
text(9,13.5, "D")
text(c(45, 135),13.5, c(50,1000))


plot(l, RuBP2[t0:t1,"All"], t = "l", ylim = c(0,100),
     ylab = expression(RuBP~(mu*mol~m^{-2}~s^{-1})), 
     xlab = "Time (s)", xaxt = "n")
axis(1, seq(0,180,30))
lines(l, RuBP2[t0:t1,"control"], col = 2)
minQSS2 = min(RuBP2[t0:t1,"QSS2"])
maxQSS2 = max(RuBP2[t0:t1,"QSS2"])
QSStime = c(0,0,90,90,180)
lines(QSStime, c(maxQSS2, minQSS2, minQSS2, maxQSS2, maxQSS2), col = 3)
abline(v =(which(abs(diff(PARdata[t0:t1])) > 100) - 1), lty = 2)
text(9,95, "E")
text(c(45, 135),95, c(50,1000))

plot(l, fRB2[t0:t1,"All"], t = "l", ylim = c(0,1),
     ylab = expression(fRB~(mu*mol~m^{-2}~s^{-1})), 
     xlab = "Time (s)", xaxt = "n")
axis(1, seq(0,180,30))
lines(l, fRB2[t0:t1,"control"], col = 2)
minQSS2 = min(fRB2[t0:t1,"QSS2"])
maxQSS2 = max(fRB2[t0:t1,"QSS2"])
QSStime = c(0,0,90,90,180)
lines(QSStime, c(maxQSS2, minQSS2, minQSS2, maxQSS2, maxQSS2), col = 3)
abline(v = (which(abs(diff(PARdata[t0:t1])) > 100) - 1), lty = 2)
text(9,0.95, "F")
text(c(45, 135),0.95, c(50,1000))

dev.off()



















# Effect of frequency -----------------------------------------------------

# Simulate control model at different frequencies
# Calculata accumulated photosynthesis

lightflecksSumA = function(model, period, PARs) {
  
  PAR1 = PARs[1]
  PAR2 = PARs[2]
  
  # CVODE settings
  model$set_settings(c("atol","rtol","maxsteps","maxerr","maxnonlin","maxconvfail","minimum"), 
                     c(1e-10,1e-8,5e4,20,20,20, -1e-6))
  model$set_settings(c("silent","positive", "force_positive"), c(TRUE, TRUE,TRUE))
  model$set_settings("maxtime", 1000)
  
  # Assign parameters from literature to the model
  filter = which(names(parameters) %in% names(model$Parameters$Values))
  model$set_parameters(names(parameters[filter]),unname(parameters[filter])) 
  
  # Simulate a LiCOR with 10% red - 90% blue
  model$set_forcings("Ib", cbind(c(0,1), c(PAR1,PAR1)*0.1))
  model$set_forcings("Ir", cbind(c(0,1), c(PAR1,PAR1)*0.9))
  model$set_forcings("Ig", cbind(c(0,1), c(0,0)))
  model$set_forcings("CO2R", cbind(c(0,1), c(400,400)))
  model$set_forcings("H2OR", cbind(c(0,1), c(20,20)))
  model$set_forcings("Ta", cbind(c(0,1), c(298.15, 298.15)))
  model$set_forcings("Tl", cbind(c(0,1), c(298.15, 298.15)))
  model$set_states("Ci", 400)
  model$set_states("Cc", 400)
  model$set_states("Ccyt", 400)
  tryCatch(model$set_states("PR", 25), error = function(x) NULL)
  
  # Calculate steady-state
  model$set_time(c(0,1800))
  steadyState = cvode(model)[2,names(model$States$Values)]
  model$set_states(names(steadyState), steadyState)
  model$set_states("sumA", 0)
  
  # Simulate Transient - Square wave light determined by period
  dt = min(1e-2, period/10)
  timeI = sort(c(seq(0,3600, by = period), seq(dt,3600 + dt, by = period)))
  model$set_forcings("Ib", cbind(timeI, c(rep(c(PAR1,PAR2,PAR2,PAR1)*0.1, times = 3600/period/2), PAR1*0.1,PAR2*0.1)))
  model$set_forcings("Ir", cbind(timeI,  c(rep(c(PAR1,PAR2,PAR2,PAR1)*0.9, times = 3600/period/2), PAR1*0.9,PAR2*0.9)))
  model$set_time(seq(0,3600,by = 1))
  diff(cvode(model)[c(1800,3600),"sumA"])/1800
}



if(.Platform$OS.type == "windows") {
  cl <- makeCluster(8)
} else {
  cl <- makeForkCluster(8)
}
registerDoParallel(cl)

periods = c(0.1, 0.3, 0.5, 1, 2, 3,  5, 10, 15, 30, 60, 90, 150, 225, 300)
LF1000 = foreach(period = periods, .packages = "MiniModel") %dopar% {
  lightflecksSumA(generate_MiniModel_model(), period, c(50,1000))
}
LF800 = foreach(period = periods, .packages = "MiniModel") %dopar% {
  lightflecksSumA(generate_MiniModel_model(), period, c(50,800))
}
LF150_1000 = foreach(period = periods, .packages = "MiniModel") %dopar% {
  lightflecksSumA(generate_MiniModel_model(), period, c(150,800))
}
LF600 = foreach(period = periods, .packages = "MiniModel") %dopar% {
  lightflecksSumA(generate_MiniModel_model(), period, c(50,600))
}
LF400 = foreach(period = periods, .packages = "MiniModel") %dopar% {
  lightflecksSumA(generate_MiniModel_model(), period, c(50,400))
}

stopCluster(cl)

AQSS1000 = mean(predict(APAR, newdata = data.frame(PAR = c(50,1000))))
AQSS800 = mean(predict(APAR, newdata = data.frame(PAR = c(50,800))))
AQSS600 = mean(predict(APAR, newdata = data.frame(PAR = c(50,600))))
AQSS400 = mean(predict(APAR, newdata = data.frame(PAR = c(50,400))))
# Check that average PAR is always the same
Aconst = predict(APAR, newdata = data.frame(PAR = mean(c(50,1000))))
png(file = "LightfleckDuration.png", width = 7, height = 4, units = "in",
    bg = "transparent", res = 1000)
par(mfrow = c(1,1), las = 1, yaxs = "i", xaxs = "i", mar = c(4,5.5,0.5,0.7),
    cex.axis = 1.2, cex.lab = 1.2, lwd = 1.2)
plot(periods, unlist(LF1000), ylim = c(4,13), log = "x", t = "o",
     xlab = "Lightfleck duration (s)",
     ylab = expression(bar(A)~(mu*mol~m^{-2}~s^{-1})))
text(100, 5, "Flashing dynamic")
abline(h = Aconst)
text(100, 12, "Constant light")
abline(h = AQSS1000, lty = 2)
text(90, 8, "Flashing steady-state")
dev.off()


png(file = "FrequencyProposal.png", width = 5, height = 4, units = "in",
    bg = "transparent", res = 1000)
par(mfrow = c(1,1), las = 1, yaxs = "i", xaxs = "i", mar = c(4,5.5,0.5,0.7),
    cex.axis = 1.2, cex.lab = 1.2, lwd = 1.2)
plot(periods, unlist(LF150_1000), ylim = c(8,13), log = "x", t = "o",
     xlab = "Fluctuation duration (s)",
     ylab = expression(Photosynthesis~(mu*mol~m^{-2}~s^{-1})),
     xaxt = "n")
axis(1, at = c(0.1,0.5,5,50))
text(100, 9, "Fluctuating")
abline(h = predict(APAR, newdata = data.frame(PAR = mean(c(150,1000)))))
text(100, 12.7, "Constant")
dev.off()


png(file = "LRC.png", width = 5, height = 4, units = "in",
    bg = "transparent", res = 1000)

par(mfrow = c(1,1), las = 1, yaxs = "i", xaxs = "i", mar = c(4,5.5,0.5,1.2),
    cex.axis = 1.2, cex.lab = 1.2, lwd = 1)

plot(PARdata, predict(APAR, newdata=data.frame(PAR = PARdata)), t = "l",
     ylab = expression(Photosynthesis~(mu*mol~m^{-2}~s^{-1})),
     xlab = expression(Light~intensity~(mu*mol~m^{-2}~s^{-1})),
     ylim = c(-1,13))

# Explain the non-linear effect
A100 = predict(APAR, newdata=data.frame(PAR = 100))
lines(c(100, 100), c(-10,A100), lty = 2, lwd = 1)
lines(c(0, 100), c(A100,A100), lty = 2, lwd = 1)

A600 = predict(APAR, newdata=data.frame(PAR = 600))
lines(c(600, 600), c(-10,A600), lty = 2, lwd = 1)
lines(c(0, 600), c(A600,A600), lty = 2, lwd = 1)

muA = mean(c(A100, A600))
A350 = predict(APAR, newdata=data.frame(PAR = 350))
lines(c(350, 350), c(-10,A350), lty = 2, lwd = 1)
lines(c(0, 350), c(A350,A350), lty = 2, lwd = 1)
lines(c(0,100), c(muA, muA))

dev.off()


# Calculate cutoff for each lightfleck intensity
cut1000 = approx(unlist(LF1000) - AQSS1000, periods, 0)$y
cut800 = approx(unlist(LF800) - AQSS800, periods, 0)$y
cut600 = approx(unlist(LF600) - AQSS600, periods, 0)$y
cut400 = approx(unlist(LF400) - AQSS400, periods, 0)$y

plot(c(400,600, 800, 1000), c(cut400, cut600, cut800, cut1000))


png(file = "InductionRelaxation.png", width = 6, height = 3, units = "in",
    bg = "transparent", res = 1000)
par(mar = c(0.5,0.5,0.5,0.5), bg = "transparent")
with(control, plot(time, Photo, t = "l", lwd = 4, col = "darkgreen", xaxt = "n",
                   yaxt = "n", bg = "transparent", bty = "n"))
dev.off()




