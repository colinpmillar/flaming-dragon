
### Colin Millar 25 December 2012
### based on code written by Johh Simmonds (January 2011)
### PLEASE NOTE THIS IS DEVELOPING CODE.
## R code used to evaluate - single stock MSY considerations
# the code fits three S-R relationships
# Bev-Holt, Hockey Stick, Ricker
# stock data are loaded as FLR objects
# years used to define selection, weights etc are used defined

library(colorout)
library(MASS)
library(FLCore)
source("SimulateRefPts2.R")
source("myMH.R")

############### Read in stock assess information prepared as FLR stock objects.
load("ewg1219.RData")

# do SRR fits
fit <- fitModels(SarEWG)

# check the fits
par(mfrow = c(2,2))
LLmodplot(fit)
SRplot(fit)

# also usefull are - note parameterisation may make plots look odd
plot(fit $ HS)
plot(fit $ RK)
plot(fit $ BH)


sim <- EqSim(fit1)



rfpts <- Eqplot(sim $ s, fit $ stk, sim $ Blim)
round(rfpts, 2)

save(fit, sim, file = "fullDataSims2.RData")

png("fullDataPlots2%01d.png", width = 1200, height = 1200, pointsize = 20)

for (i in 1:4) {
  sim <- get(paste0("sim", i))
  fit <- get(paste0("fit", i))

  par(mfrow = c(2,2))
  LLmodplot(fit $ res, fit $ data)
  SRplot(fit)
  Eqplot(sim $ s, fit $ stk, sim $ Blim)
}

dev.off()


