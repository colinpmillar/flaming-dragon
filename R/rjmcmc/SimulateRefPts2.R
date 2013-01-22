
##############################################
# utility functions
##############################################

loader <- function(p)
{
  if (p==0) cat("0%                       50%                     100%\n")
 str <- paste0(rep(c("\r[", "=", ">", " ", "]"), c(1, floor(p*50), 1, 50 - floor(p*50), 1)), collapse = "")  
 cat(str); flush.console()
 if (floor(p) == 1) cat("\n")
}

##############################################
# stock recruit formulations
##############################################

HSr <- function(ab, ssb) {
  log(ifelse(ssb >= ab$b, ab$a * ab$b, ab$a * ssb))
}

BHr <- function(ab, ssb) {
  log(ab$a * ssb / (ab$b + ssb))
}

RKr <- function(ab, ssb) {
  log(ab$a * ssb * exp(-ab$b * ssb))
}


##############################################
# Main processing subroutines: 
#   1) fitModels
#   2) EqSim
##############################################


fitModels <- function(data, runid, delta = 1.3, Nburn = 10000, use = NULL)
{

#--------------------------------------------------------
# Fit models
#--------------------------------------------------------

  delta <- rep(delta, length = 3)
  # fit stock recruit relationships using Metropolis hasings MCMC algorithm
  RK <- BHMH(Nburn + 5000, Nburn, data, delta = delta[1], model = "ricker")
  BH <- BHMH(Nburn + 5000, Nburn, data, delta = delta[2], model = "bevholt")
  HS <- BHMH(Nburn + 5000, Nburn, data, delta = delta[3], model = "segreg")

  # transform to johns parameterisations
  BH $ b <- 1/BH $ b
  BH $ a <- BH $ a * BH $ b

  a <- HS $ b
  HS $ b <- 1/(HS $ a * HS $ b)
  HS $ a <- a

#--------------------------------------------------------
# get posterior distribution of estimated recruitment
#--------------------------------------------------------
  llik.mat <- cbind(HS = HS $ llik, RK = RK $ llik, BH = BH $ llik) 
  fits <- list(
         RHS = t(sapply(1:nrow(HS), function(i) HSr(HS[i,], sort(data $ ssb)))),
         RRK = t(sapply(1:nrow(HS), function(i) RKr(RK[i,], sort(data $ ssb)))),
         RBH = t(sapply(1:nrow(HS), function(i) BHr(BH[i,], sort(data $ ssb)))))


#--------------------------------------------------------
# make stock recruit object with correct probabilities of each model
#--------------------------------------------------------
 
  # we can remove unlikely models here....  probably not nessisary

  # construct a markov chain for the model choice
  M <- rep(NA, nrow(HS))
  # random start point
  M[1] <- sample(1:3, 1)

  for (i in 2:length(M))  
  {
    # where to go next
    oldM <- M[i-1]
    newM <- sample(setdiff(1:3, oldM), 1)
    u <- runif(1)
    A <- min(1, exp(llik.mat[i,newM] - llik.mat[i-1,oldM]))
    if (u <= A) { # Accept the proposed move:
      M[i] <- newM
    } else {
      M[i] <- oldM
    }
  }

  #TODO now we should remove burn in... 

  #TODO improve this section!
  mod <- c("HS","RK","BH")[M]
  SR <- cbind(mod, HS)
  SR[mod=="RK",-1] <- RK[mod=="RK",]
  SR[mod=="BH",-1] <- BH[mod=="BH",]
  names(SR)[-1] <- c("a","b","sigma","llik")

  list(HS = HS, RK = RK, BH = BH, fits = fits, SR = SR, data = data, stknam = runid)
}



##### simulates the equilibrium results for a population
EqSim <- function(fit, stk, 
                  Nrun = 200, # number of years to run in total
                  wt.years = c(2007, 2011), # years sample weights, sel from
                  Fscan = seq(0, 1, len = 20)) # F values to scan over
{

  btyr1 <- wt.years[1]
  btyr2 <- wt.years[2] 
  flgsel <- 0
  flgmatwt <- 0
  keep <- min(Nrun, 50)

  SR <- fit $ SR
  data <- fit $ data

  # forecast settings (mean wt etc)
  stk.win <- window(stk, start = btyr1, end = btyr2)
  west <- stock.wt(stk.win)[drop=TRUE]
  weca <- catch.wt(stk.win)[drop=TRUE]
  sel <- harvest(stk.win)[drop=TRUE]
  Fbar <- fbar(stk.win)[drop=TRUE]
  sel <- sweep(sel, 2, Fbar, "/")

  if (flgsel == 0) { # take means of selection
    sel[] <- apply(sel, 1, mean)
  }
  if (flgmatwt==0){ # take means of wts
    west[] <- apply(west, 1, mean)
    weca[] <- apply(weca, 1, mean)
  } 

  Mat <- apply(mat(stk.win), 1, mean)[drop=TRUE]
  M <- mean(m(stk.win))
  Fprop <- mean(harvest.spwn(stk.win))
  Mprop <- mean(m.spwn(stk.win))

  # get ready for the simulations
  Nmod <- nrow(SR)
  NF <- length(Fscan)
  ages <- dims(stk) $ age

  ssby <- array(0, c(Nrun,Nmod))
  Ny <- Fy <- WSy <- WCy <- Cy <- Wy <- array(0, c(ages, Nrun, Nmod))
  rsam <- array(sample(1:ncol(weca), Nrun * Nmod, TRUE), c(Nrun, Nmod)) 
  Wy[] <- c(weca[, c(rsam)])

  # initial recruitment
  R <- mean( data $ rec)
  ssbs <- cats <- recs <- array(0, c(7, NF))
  pssb1 <- pssb2 <- array(0, NF)
  ssbsa <- catsa <- recsa <- array(0, c(NF, keep, Nmod))
  begin <- Nrun - keep + 1

  loader(0)
  for (i in 1:NF) {

    # The F value to test
    Fbar <- Fscan[i]

    # the selection patterns for the first year
    Zpre <- (Fbar * Fprop * sel[,rsam[1,]] + M * Mprop)
    Zpos <- (Fbar * (1-Fprop) * sel[,rsam[1,]] + M * (1-Mprop))
    # run Z out to age 50 ...
    Zcum <- c(0, cumsum(Fbar * sel[c(1:ages, rep(ages, 49 - ages)), rsam[1,]] + M))
    N1 <- R * exp(- unname(Zcum))

    # set up age structure in first year for all simulations
    Ny[,1,] <- c(N1[1:(ages-1)], sum(N1[ages:50]))

    # calculate ssb in first year using a different stock.wt for each simulation
    ssby[1,] <- colSums(Mat * Ny[,1,] * west[,rsam[1,]] / exp(Zpre))

    # loop over years
    for (j in 2:Nrun) {
      # get ssb from previous year
      SSB <- ssby[j-1,]

      # predict recruitment using various models
      Ny[1,j,] <-
          exp( ifelse(SR $ mod == "HS", HSr(SR, SSB),
                 ifelse(SR $ mod == "RK", RKr(SR, SSB), 
                   BHr(SR, SSB))) + rnorm(Nmod, 0, SR $ sigma))

      # get a selection pattern for each simulation and apply this to get N
      Zpre <- Fbar * Fprop * sel[, rsam[j,]] + M * Mprop
      Fy[    , j-1, ] <- Fbar * sel[, rsam[j-1,]]
      Ny[  -1,   j, ] <- Ny[1:(ages-1), j-1, ] * exp(-Fy[1:(ages-1), j-1, ] - M)
      Ny[ages,   j, ] <- Ny[ages, j, ] + Ny[ages, j-1, ] * exp(-Fy[ages, j-1, ] - M)
      # calculate ssb and catch.n 
      ssby[j, ] <- apply(array(Mat * Ny[,j,] * west[, rsam[j,]] / exp(Zpre), c(ages, Nmod)), 2, sum)
      Cy[, j, ] <- Ny[, j-1, ] * Fy[, j-1, ] / (Fy[, j-1, ] + M) * (1 - exp(-Fy[, j-1, ] - M))
    }

    # convert to catch weight
    Cw <- Cy * Wy
    Cat <- apply(Cw, 2:3, sum)

    # summarise everything and spit out!
    quants <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    ssbs[, i]   <- quantile(ssby[begin:Nrun, ], quants)
    cats[, i]   <- quantile(Cat[begin:Nrun, ], quants)
    recs[, i]   <- quantile(Ny[1, begin:Nrun, ], quants)
    ssbsa[i, , ] <- ssby[begin:Nrun, ]
    catsa[i, , ] <- Cat[begin:Nrun, ]
    recsa[i, , ] <- Ny[1, begin:Nrun, ]

    loader(i/NF)
  }

  list(ssbs = ssbs, cats = cats, recs = recs, ssbsa = ssbsa, catsa = catsa, recsa = recsa,
       Mat = Mat, M = M, Fprop = Fprop, Mprop = Mprop, west = west, weca = weca, sel = sel,
       Fscan = Fscan)
}


