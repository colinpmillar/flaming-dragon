
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
source("plotting.R")
source("myMH.R")

############### Read in stock assess information prepared as FLR stock objects.
load("ewg1219.RData")

# choose a stock
stk <- SarEWG

# make a directory to save results
odir <- gsub(" ","-",name(stk))
# or choose own output directory ...
res <- dir.create(odir)
if (res == FALSE) cat("  Make sure you want to write over existing files.\n  Better to empty folder before proceeding.\n")

# do SRR fits
#----------------------------

data <- data.frame(ssb = c(ssb(stk))[-dim(stock.n(stk))[2]],
                   rec = c(rec(stk))[-1])
data <- data[-9,]
runid <- name(stk)

fit <- fitModels(data, runid)
save(fit, file = paste(odir, "fit.rData", sep = "/"))

# check the fits
#----------------------------
par(mfrow = c(2,2))
LLmodplot(fit)
SRplot(fit)

# run the simulations
#----------------------------
sim <- EqSim(fit, SarEWG)
save(sim, file = paste(odir, "sim.rData", sep="/"))

# get ref points and do plots
#----------------------------
load(paste(odir, "fit.rData", sep="/"))
load(paste(odir, "sim.rData", sep="/"))

rfpts <- Eqplot(sim, stk, fit, Blim = 65000)
round(rfpts, 2)








#############################################################
#TODO do for a number of scenarios
#############################################################

png("fullDataPlots2%01d.png", width = 1200, height = 1200, pointsize = 20)

for (i in 1:4) {
  sim <- get(paste0("sim", i))
  fit <- get(paste0("fit", i))

  par(mfrow = c(2,2))
  LLmodplot(fit)
  SRplot(fit)
  Eqplot(sim $ s, fit $ stk, sim $ Blim)
}

dev.off()


