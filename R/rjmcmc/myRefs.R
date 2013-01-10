
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
source("SimulateRefPts.R")
source("myMH.R")

############### Read in stock assess information prepared as FLR stock objects.
load("ewg1219.RData")

# do SRR fits
fit1 <- fitModels(SarEWG)
fit2 <- fitModels(SarEWG.GFCM)
fit3 <- fitModels(SarGFCM)
fit4 <- fitModels(AncEWG, Nburn = 20000, delta = 1.4)

# check the fits
fit <- fit1
par(mfrow = c(2,2))
LLmodplot(fit $ res, fit $ data)
SRplot(fit)

# also usefull are - note parameterisation may make plots look odd
plot(fit $ HS)
plot(fit $ RK)
plot(fit $ BH)


sim1 <- EqSim(fit1)
sim2 <- EqSim(fit2)
sim3 <- EqSim(fit3)
sim4 <- EqSim(fit4)

sim <- sim4
fit <- fit4
rfpts <- Eqplot(sim $ s, fit $ stk, sim $ Blim)
round(rfpts, 2)

save(fit1, fit2, fit3, fit4, sim1, sim2, sim3, sim4, file = "fullDataSims.RData")

png("fullDataPlots%01d.png", width = 1200, height = 1200, pointsize = 20)

for (i in 1:4) {
  sim <- get(paste0("sim", i))
  fit <- get(paste0("fit", i))

  par(mfrow = c(2,2))
  LLmodplot(fit $ res, fit $ data)
  SRplot(fit)
  Eqplot(sim $ s, fit $ stk, sim $ Blim)
}

dev.off()



############### remove some recruits ...

# remove higest recruitment
SarEWGa <- SarEWG
stock.n(SarEWGa)[1,which.max(rec(SarEWGa))] <- NA

# remove higest recruitment
SarEWG.GFCMa <- SarEWG.GFCM
stock.n(SarEWG.GFCMa)[1,which.max(rec(SarEWG.GFCMa))] <- NA

# remove two highest recruitments
AncEWGa <- AncEWG
stock.n(AncEWGa)[1,order(rec(AncEWG), decreasing = TRUE)[1:2]] <- NA


# do SRR fits
fit1a <- fitModels(SarEWGa)
fit2a <- fitModels(SarEWG.GFCMa)
fit4a <- fitModels(AncEWGa)

sim1a <- EqSim(fit1a)
sim2a <- EqSim(fit2a)
sim4a <- EqSim(fit4a)

save(fit1a, fit2a, fit4a, sim1a, sim2a, sim4a, file = "clippedDataSims.RData")

png("clippedDataPlots%01d.png", width = 1200, height = 1200, pointsize = 20)

for (i in c(1,2,4)) {
  sim <- get(paste0("sim", i, "a"))
  fit <- get(paste0("fit", i, "a"))

  par(mfrow = c(2,2))
  LLmodplot(fit $ res, fit $ data)
  SRplot(fit)
  Eqplot(sim $ s, fit $ stk, sim $ Blim)
}

dev.off()



