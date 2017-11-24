
# Load libraries ----------------------------------------------------------
library(MiniModel)
library(ggplot2)
library(readr)
library(dplyr)
library(doParallel)
library(Hmisc)



# Irradiance transients - CO2 assimilation and gs -------------------------


# Retrieve the induction data
load("Intermediate/transients.RData")
data = filter(transients_all, Genotype == "col", TransientType %in% paste0("transient", 2:3), ID == "h")




# Figure
png("Output/figureMeasuredTransients.png", width = 14, height = 14, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(2,2), mar = rep(0,4), oma = c(4,4.5,0.5,0.6), yaxs = "i", xaxs = "i", las = 1)

with(filter(data, TransientType == "transient2", Obs > 61), {
  plot((Time - Time[1])/60, Photo, t = "l", xaxt = "n", ylim = c(0,15), xlim = c(0,35), yaxt = "n")
  axis(2, seq(0,15,2.5))
  mtext(expression(A~(mu*mol~m^{-2}~s^{-1})), side = 2, line = 2.8, las = 3, cex = 0.90)
})
text(2,14,"A")

with(filter(data, TransientType == "transient3", Obs > 61), {
  plot((Time - Time[1])/60, Photo, t = "l", yaxt = "n", xaxt = "n", ylim = c(0,15), xlim = c(0,25))
})
text(2*25/35,14,"B")

with(filter(data, TransientType == "transient2", Obs > 61), {
  plot((Time - Time[1])/60, oCond, t = "l", ylim = c(0,0.3), xlim = c(0,35), yaxt = "n")
  axis(2, seq(0, 0.25, 0.05))
  mtext(expression(g[s]~(mol~m^{-2}~s^{-1})), side = 2, line = 2.8, las = 3, cex = 0.90)
  mtext("Time (min)", side = 1, line = 2.5, las = 1, cex = 0.90, outer = TRUE)
})
text(2,0.28,"C")

with(filter(data, TransientType == "transient3", Obs > 61), {
  plot((Time - Time[1])/60, oCond, t = "l", yaxt = "n", ylim = c(0,0.3), xlim = c(0,25), xaxt = "n")
  axis(1, seq(5,25,5))
})
text(2*25/35,0.28,"D")

dev.off()



# Example lightflecks -----------------------------------------------------

# Plot A during lightfleck, plus irradiance overimposed on it

load("Intermediate/lightflecks.RData")
data = filter(lf_data_mean, Genotype == "col", Amplitude == 500, Obs < 62 + 1200)


png("Output/figureMeasuredLF.png", width = 10, height = 8, units = "cm", 
    res = 600, bg = "white", antialias = "default")

par(mfrow = c(1, 1), mar = c(4,4.0,0.5,4.0), yaxs = "i", xaxs = "i", las = 1, mgp = c(2,1,0))
with(data[-1, ], {
  plot(1, 1, ylim = c(0,13.34), xlim = c(0,21.5), t = "n", ylab = expression(A~(mu*mol~m^{-2}~s^{-1})),
       xlab = "Time (min)")
  lines((Time - Time[1])/60, PAR/45, col = 2, lwd = 0.5)
  lines((Time - Time[1])/60, Photo, lwd = 0.5)
})
axis(4, seq(0,600/45,100/45), labels = c("0", "100", "200", "300", "400", "500", "600"), col = 2, col.axis = 2)
mtext(expression(I~(mu*mol~m^{-2}~s^{-1})), side = 4, line = 2.8, las = 3, col = 2)

dev.off()