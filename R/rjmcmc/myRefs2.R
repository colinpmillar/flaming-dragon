
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


doOne <- function(runid)
{
# SGMED best assessment
if (runid == "Sardine SGMED") {
  stk <- SarEWG.GFCM
  range(stk)[c("minfbar","maxfbar")] <- c(1,4)
  data <- data.frame(ssb  = c(ssb(stk))[-dim(stock.n(stk))[2]],
                     rec  = c(rec(stk))[-1],
                     year = dimnames(catch(stk)) $ year[-1])
}

# SGMED best assessment with 1984 recruitment removed
if (runid == "Sardine SGMED outlier removed") {

  stk <- SarEWG.GFCM
  range(stk)[c("minfbar","maxfbar")] <- c(1,4)
  data <- data.frame(ssb  = c(ssb(stk))[-dim(stock.n(stk))[2]],
                     rec  = c(rec(stk))[-1],
                     year = dimnames(catch(stk)) $ year[-1])

  data <- data[-9,]
}

# GFCM best assessment
if (runid == "Sardine GFCM") {

  stk <- SarGFCM
  range(stk)[c("minfbar","maxfbar")] <- c(1,4)
  data <- data.frame(ssb  = c(ssb(stk))[-dim(stock.n(stk))[2]],
                     rec  = c(rec(stk))[-1],
                     year = dimnames(catch(stk)) $ year[-1])
}

# Full ICA assessment
if (runid == "Sardine Long ICA") {

  stk <- SarEWG
  range(stk)[c("minfbar","maxfbar")] <- c(1,4)
  data <- data.frame(ssb  = c(ssb(stk))[-dim(stock.n(stk))[2]],
                     rec  = c(rec(stk))[-1],
                     year = dimnames(catch(stk)) $ year[-1])
}

# SGMED anchovy assessment
if (runid == "Anchovy SGMED") {
  stk <- AncEWG
  range(stk)[c("minfbar","maxfbar")] <- c(1,3)
  data <- data.frame(ssb  = c(ssb(stk)),
                     rec  = c(stock.n(stk)["0"]),
                     year = dimnames(catch(stk)) $ year)
}

# SGMED assessment with
if (runid == "Anchovy SGMED outlier removed") {

  stk <- AncEWG
  range(stk)[c("minfbar","maxfbar")] <- c(1,3)
  data <- data.frame(ssb  = c(ssb(stk)),
                     rec  = c(stock.n(stk)["0"]),
                     year = dimnames(catch(stk)) $ year)

  data <- subset(data, ssb < 550000)
}

if (runid == "Anchovy SGMED age 0 removed") {
  stk <- AncEWG[paste(1:5),]
  range(stk)[c("minfbar","maxfbar")] <- c(1,3)
  data <- data.frame(ssb  = c(ssb(stk))[-dim(stock.n(stk))[2]],
                     rec  = c(stock.n(stk)["1"])[-1],
                     year = dimnames(catch(stk)) $ year[-1])
}


# save data
#----------------------------
odir <- gsub(" ","-",runid)

# or choose own output directory ...
res <- dir.create(odir)
if (res == FALSE) cat("  Make sure you want to write over existing files.\n  Better to empty folder before proceeding.\n")

save(odir, runid, data, stk, file = paste(odir, "data.rData", sep = "/"))

                  
# do SRR fits
#----------------------------
fit <- fitModels(data, runid)
save(fit, file = paste(odir, "fit.rData", sep = "/"))

# run the simulations
#----------------------------
sim <- EqSim(fit, stk, Fscan = seq(0, 2, length = 40))
save(sim, file = paste(odir, "sim.rData", sep="/"))


# plot the data
#----------------------------
png(paste0(odir, "/data.png"), 800, 800)                 
with(data,{
  ssb.scl <- 1000
  rec.scl <- 100000
  plot(ssb/ssb.scl, rec/rec.scl, type = "n", las = 1, 
       xlab  = sprintf("SSB (%is)", ssb.scl), 
       ylab =  sprintf("Recruits (%is)", rec.scl),
       main = paste(runid, "stock-recruit data"))
  lines(ssb/ssb.scl, rec/rec.scl, type = "b", cex = 3)
  text(x = ssb/ssb.scl, y = rec/rec.scl, label = substring(year, 3, 4), cex = 0.7)

})
dev.off()
   
# plot the fits
#----------------------------

png(paste0(odir, "/HS-fit.png"), 800, 800)
plot(fit $ HS)
dev.off()

png(paste0(odir, "/RK-fit.png"), 800, 800)
plot(fit $ RK)
dev.off()

png(paste0(odir, "/BH-fit.png"), 800, 800)
plot(fit $ BH)
dev.off()

png(paste0(odir, "/model-fits.png"), 800, 800)
par(mfrow = c(2,2))
LLmodplot(fit)
SRplot(fit)
dev.off()


# get ref points and do plots
#----------------------------
load(paste(odir, "fit.rData", sep="/"))
load(paste(odir, "sim.rData", sep="/"))

png(paste0(odir, "/refpts.png"), 800, 800)
rfpts <- Eqplot(sim, stk, fit, Blim = 0.3 * max(fit $ data $ ssb))
dev.off()

write.table(as.data.frame(t(c(round(rfpts[1:2]), round(rfpts, 2)[-c(1:2)]))), file = paste0(odir, "/refpts.csv"), 
            sep = ",", row.names = FALSE)

}


# choose a scenario

runids <- c("Sardine SGMED", "Sardine SGMED outlier removed", "Sardine GFCM", "Sardine Long ICA", 
            "Anchovy SGMED", "Anchovy SGMED outlier removed", "Anchovy SGMED age 0 removed")

library(multicore)
tmp <- mclapply(runids, doOne)

tmp <- mclapply(runids, doOne)


