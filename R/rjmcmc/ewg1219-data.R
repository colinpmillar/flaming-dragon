

library(FLCore)

getStk <- function(dir) {
  stk <- readFLStock(paste0(dir, "7511index.dat"))
  catch.n(stk) <- readVPAFile(paste0(dir, "7511cn.dat"))
  catch.wt(stk) <- readVPAFile(paste0(dir, "7511cw.dat"))
  catch(stk) <- computeCatch(stk)

  stk
}

addICA <- function(dir, stk) {

  file <- paste0(dir,"ica.vie")
  ica.vie <- readLines(file)
  floc <- grep("FISHING MORTALITY AT AGE AND YEAR", ica.vie)
  nloc <- grep("NUMBERS AT AGE AND YEAR", ica.vie)
  yrs <- scan(file, skip = 1, n = 2)

  harvest(stk)[, paste(yrs[1]:yrs[2])] <- c(t(read.table(file, skip = floc, nrows = diff(yrs)+1)))
  stock.n(stk)[, paste(yrs[1]:yrs[2])] <- c(t(read.table(file, skip = nloc, nrows = diff(yrs)+1)))
  stock(stk) <- computeStock(stk)
  
  stk
}


Sdir <- "ewg-results/Sardine/S"
Adir <- "ewg-results/Anchovy/A"

SLdir <- "ewg-results/Sardine/LongRun/"
SSdir <- "ewg-results/Sardine/ShortRun/"
ALdir <- "ewg-results/Anchovy/"


# get stocks for model fitting
SarEWG <- addICA(SLdir, getStk(Sdir))
name(SarEWG) <- "Sardine ewg"
SarEWG.GFCM <- addICA(SSdir, SarEWG)
name(SarEWG.GFCM) <- "Sardine ewg + gfcm"

SarGFCM <- window(SarEWG.GFCM, start = 2000, end = 2011)
name(SarGFCM) <- "Sardine gfcm"

#plot(FLStocks(long = SarEWG, short = SarEWG.GFCM))


AncEWG <- addICA(ALdir, getStk(Adir))
name(AncEWG) <- "Anchovy ewg"

plot(AncEWG)

# save!
save(SarEWG, SarGFCM, SarEWG.GFCM, AncEWG, file = "ewg1219.RData")


