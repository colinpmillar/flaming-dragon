


fitModels <- function(stk, delta = 1.3, Nburn = 10000)
{

  data <- data.frame(ssb = c(ssb(stk))[-dim(stock.n(stk))[2]],
                   rec = c(rec(stk))[-1])

  # remove any NAs
  data <- data[apply(data, 1, function(x) all(!is.na(x))), ]


#########################################################
# Fit models
#######################################################

  # fit stock recruit relationships using Metropolis hasings MCMC algorithm
  RK <- BHMH(Nburn + 5000, 5000, data, delta = delta, model = "ricker")
  BH <- BHMH(Nburn + 5000, 5000, data, delta = delta, model = "bevholt")
  HS <- BHMH(Nburn + 5000, 5000, data, delta = delta, model = "segreg")

  # transform to johns parameterisations
  BH $ b <- 1/BH $ b
  BH $ a <- BH $ a * BH $ b

  a <- HS $ b
  HS $ b <- 1/(HS $ a * HS $ b)
  HS $ a <- a

#########################################################
# get posterior distribution of estimated recruitment
#######################################################
  res <- 
   list(llHS = HS $ llik, llRK = RK $ llik, llBH = BH $ llik, 
         RHS = t(sapply(1:nrow(HS), function(i) HSr(HS[i,], sort(data $ ssb)))),
         RRK = t(sapply(1:nrow(HS), function(i) RKr(RK[i,], sort(data $ ssb)))),
         RBH = t(sapply(1:nrow(HS), function(i) BHr(BH[i,], sort(data $ ssb)))),
         Blim = median(HS $ b), PBlim = hist(HS $ b, 30, plot=FALSE),
         Stknam = name(stk))


#########################################################
# make stock recruit object with correct probabilities of each model
#######################################################
  SR <- SRobject(getProb(res, 1000), 
                 HS[1:1000,c("a","b","cv")], 
                 RK[1:1000,c("a","b","cv")], 
                 BH[1:1000,c("a","b","cv")])

  list(HS = HS, RK = RK, BH = BH, res = res, SR = SR, stk = stk, data = data)
}



EqSim <- function(fit, Nrun = 200, # number of years to run in total 
                  wt.years = c(2006, 2011), # years sample weights, sel from
                  Fscan = seq(0.01, 1, len = 20), # F values to scan over
                  Blim = max(fit $ data $ ssb) * 0.3) # Blim - deafults to 30% of maximum Rec
{

  out <- Eqsym(fit, wt.years[1], wt.years[2], Fscan, Nrun, Blim, keep = 100) 

  list(simlist = out, wt.years = wt.years, Blim = Blim)
}







#############################################
############################################
#Main processing subroutines
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


# get posterior model probabilities
getProb <- function(res, n) {

  llHS = sort(res[[1]])[1:n]
  llRK = sort(res[[2]])[1:n]
  llBH = sort(res[[3]])[1:n]
 
  # scale likelihoods to avoid numerical overflow
  sfactor <- mean(c(llHS, llRK, llBH))

  P1t = 1/mean(exp(-llHS + sfactor)) 
  P3t = 1/mean(exp(-llBH + sfactor)) 
  P7t = 1/mean(exp(-llRK + sfactor)) 

  # HS RK BH
  c(P1t,P3t,P7t)/sum(P1t,P3t,P7t)
}


######### plot the fitted model curves

LLplot <- function(res) 
{
  P1a = sort(res[[1]])
  P3a = sort(res[[2]])
  P7a = sort(res[[3]])
  plot(P1a, type = "n", 
       main = paste("LL of Bayes model:", res $ Stknam), 
       xlab = "model order", ylab = "log likelihood", 
       ylim = range(P1a,P3a,P7a))
  lines(P1a, lty = 1, col = 1)
  lines(P3a, lty = 2, col = 2)
  lines(P7a, lty = 3, col = 3)
  legend(x = "bottomright", legend = c("HS","RK","BH"), 
         lty = c(1,2,3), col = c(1,2,3)) 
}

LLmodplot <- function (res, data)
{
  Rec <- data $ rec[order(data $ ssb)]
  SSB <- sort(data $ ssb)

  for (mod in 1:3) 
  { 
    model <- c("Hockey Stick", "Ricker", "Beverton Holt")[mod]
    Rsym <- exp(res[[3 + mod]])

    plot(SSB, Rec, xlab = 'SSB', ylab = 'Recruits', 
         main = paste(res $ Stknam, model), 
         ylim = c(0, max(Rec)), xlim = c(0, max(SSB)))
    for (i in sample(1:nrow(res[[4]]), 1000)) {
      lines(SSB, Rsym[i,],col = paste0(grey(0), "10") )
    }
    for (i in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
      lines(SSB, apply(Rsym, 2, quantile, i), col = 2)
    }
    lines(SSB, Rsym[which.max(res[[mod]]),], col=1, lwd=2)
    points(SSB, Rec, pch = 19, col = "darkblue")  
  }
}


#### subroutine to creat modeset of 1000 models with set defined in 'mod'.
SRobject <- function(ps, HS, RK, BH) 
{
  mod <- sample(c("HSL","RKL","BHL"), 1000, replace = TRUE, prob = ps)
  pp <- 1:1000
  RKLp <- pp[mod=="RKL"]
  BHLp <- pp[mod=="BHL"]

  ###### without truncation
  Modset <- cbind(mod, HS)
  Modset[RKLp,2:4] <- RK[RKLp,1:3]
  Modset[BHLp,2:4] <- BH[BHLp,1:3]
  names(Modset)[2:4] <- c("A","B","sigma")

  Modset 
}

#######################################################
# plot out recruits without truncation  including ----- Blim estimates
SRplot <- function (fit) 
{

  Modset <- fit $ SR
  data <- fit $ data
  res <- fit $ res

  SSBO <- data $ ssb
  RecO <- data $ rec

  StkNam <- res $ Stknam
  Blim <- res $ Blim
  PBlim <- res $ PBlim


  maint=paste(StkNam,"Sym SR")
  mn=length(SSBO)
  SSBs=seq(1:102000)
  Recs=seq(1:102000)
  minSSB=min(c(SSBO,max(SSBO)*0.05))
  maxSSB=max(SSBO)*1.1
  maxrec=max(RecO*1.5)
  plot(SSBO,RecO,xlim=c(0,maxSSB),ylim=c(0,maxrec),type="p",pch=19,col=10,xlab="SSB ('000 t)",ylab="Recruits",main=StkNam)
  for (i in 1:1000) {
    SSB = runif(1020, minSSB, maxSSB)
    #SSB = rep(SSBO,40) ## to match SSB values to original
    if (Modset $ mod[i] =="HSL") {
      mu=log(Modset$A[i]*Modset$B[i] * (SSB>=Modset$B[i]) + Modset$A[i]*SSB*(SSB<Modset$B[i]))
      R2=rnorm(1020,0,Modset$sigma[i])
      Rec=exp(mu+R2[1:1020])
    } else
    if (Modset $ mod[i] == "RKL") {
      mu=log(Modset$A[i]*SSB*exp(-Modset$B[i]*SSB))
      R2=rnorm(1020,0,Modset$sigma[i])
      Rec=exp(mu+R2[1:1020])
    } else
    if (Modset $ mod[i] == "BHL") {
      mu <-log(Modset$A[i]*SSB/(Modset$B[i]+SSB))
      R2=rnorm(1020,0,Modset$sigma[i])
      Rec=exp(mu+R2[1:1020])
    }
    points(SSB[1:mn],Rec[1:mn],type="p",pch=20,col=1,cex=0.0625)
    SSBs[((i-1)*1020+1):(i*1020)]=SSB
    Recs[((i-1)*1020+1):(i*1020)]=Rec
  }
  points(SSBO,RecO,type="p",pch=19,col=10,cex=1.25)

  step=maxSSB*0.05/1.1
 
  up=seq(minSSB,maxSSB,step)
  lw=seq(minSSB,maxSSB,step)
  md=seq(minSSB,maxSSB,step)
  ssb=seq(minSSB+step/2,maxSSB-step/2,step)
  for (j in 1:(length(up)-1)){
    up[j]=quantile(Recs[((SSBs>up[j])*(SSBs<up[j+1]))>0],probs=.95,na.rm=TRUE)
    lw[j]=quantile(Recs[((SSBs>lw[j])&(SSBs<lw[j+1]))>0],probs=.05,na.rm=TRUE)
    md[j]=quantile(Recs[((SSBs>md[j])&(SSBs<md[j+1]))>0],probs=.5,na.rm=TRUE)
  }
  lines(ssb,md[1:(length(up)-1)],col=7,lwd=3)
  lines(ssb,up[1:(length(up)-1)],col=4,lwd=3)
  lines(ssb,lw[1:(length(up)-1)],col=4,lwd=3)

  # plot Blim
  #lines(PBlim$mids,0.1*maxrec/max(PBlim$counts)*PBlim$counts,col=5,lwd=2)
  lines(c(Blim,Blim),c(0,maxrec),col=5,lwd=2)
}




loader <- function(p)
{
  if (p==0) cat("0%                       50%                     100%\n")
 str <- paste0(rep(c("\r[", "=", ">", " ", "]"), c(1, floor(p*50), 1, 50 - floor(p*50), 1)), collapse = "")  
 cat(str); flush.console()
 if (floor(p) == 1) cat("\n")
}

##### simulates the equilibrium results for a population
Eqsym <- function (fit, btyr1, btyr2, 
                   Fscan, Nrun, Blim, Bpa = Blim * 1.4,  
                   flgsel = 1, flgmatwt = 1, keep = min(Nrun, 50)) 
{

###########################

  SRO <- fit $ SR
  stk <- fit $ stk
  data <- fit $ data

  # forecast settings (mean wt etc)
  stk.win <- window(stk, start = btyr1, end = btyr2)
  west <- stock.wt(stk.win)[drop=TRUE]
  weca <- catch.wt(stk.win)[drop=TRUE]
  sel <- harvest(stk.win)[drop=TRUE]
  Fbar <- fbar(stk.win)[drop=TRUE]
  sel <- sweep(sel, 1, Fbar, "/")

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
  Nmod <- dim(SRO)[1]
  NF <- length(Fscan)
  ages <- dims(stk) $ age

  ssby <- array(0, c(Nrun,Nmod))
  Ny <- Fy <- WSy <- WCy <- Cy <- Wy <- array(0, c(ages, Nrun, Nmod))
  rsam <- array(sample(1:length(Mat), Nrun * Nmod, TRUE), c(Nrun, Nmod)) 
  Wy[] <- c(weca[, rsam])

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
    # run Z out to 50 ...
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
        with(SRO,
          ifelse(mod == "HSL",
            A * B * (SSB >= B) + A * SSB * (SSB < B),
          ifelse(mod == "RKL",
            A * SSB * exp(-B * SSB),
            A * SSB / (B + SSB)
          )) * exp(rnorm(Nmod, 0, sigma)))

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
    ssbs[,i]   <- quantile( ssby[begin:Nrun,], quants)
    cats[,i]   <- quantile(Cat[begin:Nrun,], quants)
    recs[,i]   <- quantile(Ny[1,begin:Nrun,], quants)
    pssb1[i]   <- sum( ssby[begin:Nrun,] > Blim) / Nmod  /(Nrun-begin+1)
    pssb2[i]   <- sum( ssby[begin:Nrun,] > Bpa) / Nmod / (Nrun-begin+1)
    ssbsa[i,,] <- ssby[begin:Nrun, ]
    catsa[i,,] <- Cat[begin:Nrun, ]
    recsa[i,,] <- Ny[1, begin:Nrun, ]

    loader(i/NF)
  }  # new Fbar

  catm <- apply(catsa, 1, mean)
  maxcatm <- which.max(catm)
  catsam <- apply(catsa,c(1,3),mean)
  maxpf <- apply(catsam,2,which.max)
  fmsy <- Fscan[maxpf]
  v <- hist(fmsy, breaks = c(Fscan-Fscan[1]/2,max(Fscan)+Fscan[1]/2), plot = FALSE)

  pp1 <- max(which(pssb1>.95))
  grad <- (Fscan[pp1+1]-Fscan[pp1])/(pssb1[pp1+1]-pssb1[pp1])
  flim <- Fscan[pp1]+grad*(0.95-pssb1[pp1])

  if (any(cumsum(v$counts) > Nmod / 2 - 1)) {
    upp <- min(cumsum(v$counts)[ cumsum(v$counts) > Nmod / 2 - 1])
    uppp <- which(cumsum(v$counts)==upp)
  } else upp <- NA

  if (any(cumsum(v$counts) < Nmod / 2 - 1)) {
    dwn <- max(cumsum(v$counts)[ cumsum(v$counts) < Nmod / 2 - 1])
    dwnp <- which(cumsum(v$counts)==dwn)
  } else dwn <- NA

  if (!is.na(upp) & !is.na(dwn)) {
    vcum <- (dim(Cat)[2]/2-dwn)*(Fscan[uppp]-Fscan[dwnp])/(upp-dwn)+Fscan[dwnp]   # 50%  on dist msy
  } else vcum <- NA
  vmode <- which.max(v$counts)                                                    # mode of dist msy
  msym <- sum(Fscan * v$counts) / sum(v$counts)                                   # mean of dist msy

  list(Fscan,data $ rec, data $ ssb,recs,ssbs,cats,catm,flim,pssb1,pssb2,grad,
       v,vcum,vmode,msym,fmsy,maxpf,catsam,ssbsa,catsa,recsa)
} 



######## Creates equilibrium plots

Eqplot <- function (symlist, stk, Blim, Bpa = 1.4 * Blim)
{
  Nmod <- 1000 # fixed

  maint <- name(stk)  

  Fscan <- symlist[[1]]
  RecO <- symlist[[2]]
  SSBO <- symlist[[3]]

  Catchs <- catch(stk)[, 1:length(SSBO), drop = TRUE]
  FbarO <- fbar(stk)[, 1:length(SSBO), drop = TRUE]

  recs <- symlist[[4]]
  ssbs <- symlist[[5]]
  cats <- symlist[[6]]
  catm <- symlist[[7]]
  flim <- symlist[[8]]
  pssb1 <- symlist[[9]]
  pssb2 <- symlist[[10]]
  grad <- symlist[[11]]
  v <- symlist[[12]]
  vcum <- symlist[[13]]
  vmode <- symlist[[14]]
  msym <- symlist[[15]]
  fmsy <- symlist[[16]]
  ssbsa <- symlist[[19]]
  catsa <- symlist[[20]]
  recsa <- symlist[[21]]

  set <- dim(ssbsa)
  Nmod <- set[3]
  keep <- set[2]
  NF <- set[1]
  pssb3 <- apply(ssbsa>Blim,1,sum)/Nmod/keep
  pp1 <- max(which(pssb3>.50))
  grad <- (Fscan[pp1+1]-Fscan[pp1])/(pssb3[pp1+1]-pssb3[pp1])
  flim5 <- Fscan[pp1]+grad*(0.5-pssb3[pp1])

  maxcatm <- which.max(catm)

  par(mfrow=c(2,2),mar=c(2.5, 2,1.5, 1),oma=c(0,0,0,0),cex.axis=0.75, tcl=0.25,mgp=c(0,0.25,0),las=1)

  xmax <- max(Fscan)

  ymax=max(recs[7,],RecO)
  plot(Fscan,recs[7,],type="l",lty=4, ylim=c(0,ymax),xlim=c(0,xmax),ylab="",xlab="")
  title(ylab="Recruitment", xlab="F 2-6",cex.lab=0.75,line=1.2,cex.main=0.75)
  mtext(text=paste(maint," a) Recruits"),cex=0.75,side=3,line=0.5)
  lines(Fscan,recs[6,],type="l",lty=3)
  lines(Fscan,recs[5,],type="l",lty=2)
  lines(Fscan,recs[4,],type="l",lty=1)
  lines(Fscan,recs[3,],type="l",lty=2)
  lines(Fscan,recs[2,],type="l",lty=3)
  lines(Fscan,recs[1,],type="l",lty=4)
  points(FbarO,RecO,pch=21,cex=.75,bg=1)
  lines(c(flim,flim),c(0,ymax),type="l",col=3)

  ymax=max(ssbs[7,])
  plot(Fscan,ssbs[7,],type="l",lty=4,ylim=c(0,ymax),xlim=c(0,xmax),ylab="",xlab="")
  title( ylab="SSB", xlab="F 2-6",cex.lab=0.75,line=1.2,cex.main=0.75)
  mtext(text="b) Spawning Stock Biomass",cex=0.75,side=3,line=0.5)
  lines(Fscan,ssbs[6,],type="l",lty=3)
  lines(Fscan,ssbs[5,],type="l",lty=2)
  lines(Fscan,ssbs[4,],type="l",lty=1)
  lines(Fscan,ssbs[3,],type="l",lty=2)
  lines(Fscan,ssbs[2,],type="l",lty=3)
  lines(Fscan,ssbs[1,],type="l",lty=4)
  lines(c(0,xmax),c(Blim,Blim))
  text(x=0.1, y=Blim*1.1,"Blim",cex=0.7)
  points(FbarO,SSBO,pch=21,cex=.75,bg=1)
  lines(c(flim,flim),c(0,ymax),type="l",col=3)



  ymax=max(cats[7,])
  plot(Fscan,cats[7,],type="l",lty=4,ylim=c(0,ymax),xlim=c(0,max(Fscan)),ylab="",xlab="")
  title(ylab="Catch", xlab="F 2-6",cex.lab=0.75,line=1.2,cex.main=0.75)
  mtext(text="c) Catch",cex=0.75,side=3,line=0.5)
  #points(fbarsa[ssbsa<BLIM],catsa[ssbsa<BLIM],pch=21,col=6,bg=6,cex=0.125)
  #points(fbarsa[ssbsa>BLIM]+.0075,catsa[ssbsa>BLIM],pch=21,col=7,bg=7,cex=0.125)
  lines(Fscan,cats[7,],type="l",lty=4)
  lines(Fscan,cats[6,],type="l",lty=3)
  lines(Fscan,cats[5,],type="l",lty=2)
  lines(Fscan,cats[4,],type="l",lty=1)
  lines(Fscan,cats[3,],type="l",lty=2)
  lines(Fscan,cats[2,],type="l",lty=3)
  lines(Fscan,cats[1,],type="l",lty=4)
  points(FbarO,Catchs,pch=21,cex=.75,bg=1)
  lines(c(flim,flim),c(0,ymax),type="l",col=3)
  lines(Fscan,catm,type="l",lty=1,col=2)
  lines(rep(Fscan[maxcatm],2),c(0,ymax),type="l",lty=1,col=5)


  plot(Fscan,1-pssb1,type="l",lty=2,ylim=c(0,1),xlim=c(0,max(Fscan)),ylab="",xlab="")
  title(ylab="Prob MSY, SSB<Bpa or Blim", xlab="F 2-6",cex.lab=0.75,line=1.2,cex.main=0.75)
  mtext(text="d) Prob MSY and Risk to SSB",cex=0.75,side=3,line=0.5)
  lines(Fscan,1-pssb2,type="l",lty=4)
  lines(c(0,xmax),c(0.05,0.05))
  lines(c(0,xmax),c(0.1,0.1))
  text(x=0.1, y=0.075,"5%",cex=0.75)
  text(x=max(Fscan[pssb2>0.5])-.05,y=0.5,"SSB<Bpa",cex=0.75)
  text(x=max(Fscan[pssb1>0.7])+.1,y=0.3,"SSB<Blim",cex=0.75)
  lines(c(flim,flim),c(0,1),type="l",col=3)
  lines(v$mids,v$counts/Nmod,type="l",col=4)
  text(x=0.2,y=0.15,"Prob of Fmsy",cex=0.75)
  lines(rep(vcum,2),c(0,1),type="l",lty=1,col=5)

  FCrash5  = Fscan[which(cats[2,]<(max(cats[2,])/20))[1]]
  FCrash50 = Fscan[which(cats[4,]<(max(cats[4,])/20))[1]]

  out <- c(Blim,Bpa,round(vcum*100)/100,Fscan[maxcatm],flim, NA,FCrash5,FCrash50)
  names(out) <- c("Blim","Bpa","MSY:median","Maxmeanland","Flim","Flim5","FCrash5","FCrash50")

  out
}


pssb <- function (symlistf,Blim,Bpa,stknam)   {
  maintfn=paste(stknam,'ssb sum')
  maxc=max(symlistf[[7]])
  intc=which(symlistf[[7]]>0.95*maxc)
  lc=intc[1]
  mc=which.max(symlistf[[7]])
  uc=intc[length(intc)]
  npoint=length(symlistf[[19]][mc,,])
  qq= hist(symlistf[[19]][mc,,],100,plot=FALSE)
  plot(qq$breaks,c(0,cumsum(qq$counts)/npoint),xlim=c(0,1.3*symlistf[[5]][7,mc]),type='l',
       main=stknam,xlab='ssb',ylab='Cumulative Probability of SSB',lwd=2)
  qq= hist(symlistf[[19]][lc,,],100,plot=FALSE)
  lines(qq$breaks,c(0,cumsum(qq$counts)/npoint),type='l',lwd=1,lty=2)
  qq= hist(symlistf[[19]][uc,,],100,plot=FALSE)
  lines(qq$breaks,c(0,cumsum(qq$counts)/npoint),type='l',lwd=1,lty=2)

  lines(rep(Bpa,2),c(0,1),col=4,lwd=2)
  text('Bpa',x=Bpa,y=0.8)
  lines(rep(Blim,2),c(0,1),col=3,lwd=2)
  text('Blim',x=Blim,y=0.7)
  lines(rep(symlistf[[5]][4,mc],2),c(0,1),col=2,lwd=2)
  text('Median',x=symlistf[[5]][4,mc],y=0.9)
  text(paste('F=',Fscan[mc],sep=""),x=symlistf[[5]][4,mc],y=0.85)
  lines(rep(symlistf[[5]][4,uc],2),c(0,1),col=2,lwd=1,lty=2)
  lines(rep(symlistf[[5]][4,lc],2),c(0,1),col=2,lwd=1,lty=2)
  text(paste('F=',Fscan[uc],sep=""),x=symlistf[[5]][4,uc],y=0.95)
  text(paste('F=',Fscan[lc],sep=""),x=symlistf[[5]][4,lc],y=0.95)

}


 ##### set of rountines to give summary plotting of stock info
GROSEL <- function (stk,btyr1,btyr2) 
{
  fbarage=range(stk)[c('minfbar','maxfbar')]
  agrange=range(stk)[c('min','max')]
  agrange=seq(agrange[1],agrange[2])
  Stk.bwin=window(stk,start=btyr1,end=btyr2)
  mat=Stk.bwin@mat
  west=Stk.bwin@stock.wt
  weca=Stk.bwin@catch.wt
  sel=Stk.bwin@harvest
  Fbar=fbar(Stk.bwin)
  for (j in 1:length(Fbar)) {sel[,j]=sel[,j]/as.real(Fbar[,j])     }
  list(west,weca, sel,agrange,mat)
}


multline <- function (dat,age,num) {
  for (j in 1:dim(dat)[2]) {
    lines(age,dat[,j,],col=num,lty=num,type='l',lwd=2)
  }
}

singline <- function (dat,age,num) {
  lines(age,apply(dat,1,mean),col=num,lty=num,type='l',lwd=2)
}

maxminline <- function (dat,age,num) {
  lines(age,apply(dat,1,max),col=num,lty=num,type='l',lwd=2)
  lines(age,apply(dat,1,min),col=num,lty=num,type='l',lwd=2)
}



plotRefFits <- function(fit, width = 1200, height = 1200, pointsize = 30)
{

  stk <- fit $ stk
  res <- fit $ res
  SRO <- fit $ SRO

  StkNam <- name(stk)
  ###### list of S-R model names used
  ModName <- c('Hockey Stick', 'Ricker', 'Beverton Holt')

  png(paste0(gsub(" ",".",StkNam), "params%03d.png"), width = width, height = height, units = "px", pointsize = pointsize, bg = "white")

####### plot likelihood
  LLplot(res)

####### plot models
  LLmodplot(res, data)


##########################################################################################
########## use one for each model type each stock - parameter distribution plots

  Paramp(data $ ssb, data $ rec, 
         HS $ a,HS $ b, HS $ cv, which.max(res[[1]]), 
         HS $ a,HS $ b, HS $ cv, StkNam, 'Hockey Stick')  

  Paramp(data $ ssb, data $ rec, 
         RK $ a, RK $ b, RK $ cv, which.max(res[[2]]), 
         RK $ a, RK $ b, RK $ cv, StkNam,'Ricker')  

  Paramp(data $ ssb, data $ rec, 
         BH $ a, BH $ b, BH $ cv, which.max(res[[3]]), 
         BH $ a, BH $ b, BH $ cv,StkNam,'Beverton Holt')  


  dev.off()

  png(paste0(gsub(" ",".",StkNam), "srr%03d.png"), width = width, height = height, units = "px", pointsize = pointsize, bg = "white")

  #Plot simulated  S-R values to compare with observed
  SRplot(fit) 

  dev.off()


}




plotEqFits <- function(rdafile, width = 1200, height = 1200, pointsize = 30)
{

  load(rdafile)

  StkNam <- name(stk)

  png(paste0(gsub(" ",".",StkNam), "equilibrium%03d.png"), width = width, height = height, units = "px", pointsize = pointsize, bg = "white")


  pts <- Eqplot(symlistf, stk, Blim) 


########################## section with summary plots

#### summary SSB plots
  pssb(symlistf, Blim, Bpa, StkNam)

  plot(Fscan, symlistf[[7]]/max(symlistf[[7]]), lwd=2, 
       xlab='Target F', ylab = 'Rel Catch', type = 'l', main = 'Mean Catch')

  v <- symlistf[[12]]
  plot(v$mids, v$counts/max(v$counts), type='l', lwd=2,
       xlab = 'Target F', ylab = 'Rel Prob FMSY', main = 'Probability F(MSY)')

  GS = GROSEL(stk, btyr1, btyr2)

  selmax = max(GS[[3]])
  cwtmax = max(GS[[2]])
  swtmax = max(GS[[1]])
  age = c(min(GS[[4]]), max(GS[[4]]))


  plot(age, c(0,selmax), main = 'Selection relative to Fbar', 
                         xlab = 'age', ylab = 'Selection',col = 0)
  multline(GS[[3]], GS[[4]],1)


  plot(age, c(0,selmax), main = 'Selection relative to Fbar',
                         xlab = 'age', ylab = 'Mean Selection', col = 0)
  singline(GS[[3]], GS[[4]], 1)


  plot(age, c(0,selmax), main = 'Selection relative to Fbar', 
                         xlab = 'age', ylab = 'Max-Min Selection', col = 0)
  maxminline(GS[[3]], GS[[4]], 1)

  plot(age, c(0,swtmax), main = 'Weight at age in Stock',
                         xlab = 'age', ylab = 'Weight', col = 0)
  multline(GS[[1]], GS[[4]], 1)


  plot(age, c(0,cwtmax), main = 'Weight at age in Catch',
                         xlab = 'age', ylab = 'Weight', col = 0)
  multline(GS[[2]], GS[[4]], 1)


  plot(age, c(0,cwtmax), main = 'Mean Weight at age in Catch',
                         xlab = 'age', ylab = 'Weight', col = 0)
  singline(GS[[2]], GS[[4]],1)


  plot(age, c(0,cwtmax), main = 'Max-Min Weight at age in Catch',
                         xlab = 'age', ylab = 'Weight', col = 0)
  maxminline(GS[[2]], GS[[4]],1)

  plot(age, c(0,1), main = 'Mean Maturity at age', 
                    xlab = 'age', ylab = 'Fraction Mature', col = 0)
  singline(GS[[5]], GS[[4]],1)

  dev.off()

  pts
}






#### plot parameter distributions for 2 parameter models
Paramp <- function (SSB,Rec,HSa,HSb,sigmaHS,pp1,FHSa,FHSb,FHSsigma,StkName,ModName)  
{
  maint=paste(StkName," ",ModName)

  LLa=HSa[pp1]
  LLb=HSb[pp1]
  LLsigmaHS=sigmaHS[pp1]
  AA=FHSa
  BB=FHSb

 ### run this plotiing bit for each bayes model ------------ Routine borrowed and amended from Mark Payne

  AA.est=median(AA)
  BB.est=median(BB)

  # some settings normally set in the subroutine call
  n.grid=50    # this and the n below will control smoothness
  show.points=TRUE
  show.ll=TRUE
  do.contours=TRUE
  filled.contours=FALSE
  f.ages=NULL
  margin.plots=TRUE
  xlim= max(BB)
  ylim= max(AA)
  debug=FALSE
  n=40000 # set to match my full data length
  pch="."
  show.grid=TRUE
  alpha=0.05   # sets intervals on pfs
  show.estimate=TRUE
  thin=20
  TH=seq(thin/2,n,thin)

# main contouring section

  kern  <-  kde2d(BB,AA,n=n.grid)
  #Calculate cumulative distribution function
  kz    <-  as.vector(kern$z)
  ord   <-  order(kz)
  cumfrac <-  cumsum(kz[ord])/sum(kz)
  cumfrac.matrix  <-  matrix(cumfrac[rank(kz)],nrow=nrow(kern$z),ncol=ncol(kern$z))
  contour.args <- list()
  contour.args$levels <-   c(0.1,0.25,0.50,0.75,0.90)
  contour.args$labels <-   NULL
#      if(is.null(contour.args$lty))   contour.args$lty <-   c(1,1,2,3,4)
#      if(is.null(contour.args$lwd))   contour.args$lwd <-   c(1,3,1,1,1)
  if(is.null(contour.args$method))contour.args$method <-    "edge"
#    if(is.null(contour.args$labels))contour.args$labels <-    NULL
#      if(filled.contours) {      
#       do.call(filled.contour,c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,nlevels=100,color.palette=heat.colors))
#    }
  otolith.obj  <-  c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,contour.args)
# this form works but gives odd contour labeling  - alternative crude seting at end to deal with issue
# I cannot work out how to change contour.args the iff line 6 above commented out was an attemp that did not work 

  if(!show.points) { pch <- NA }
  if(margin.plots) {
    layout(matrix(c(1,4,3,2),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
    par(mar=c(0.5,0.5,0.5,0.5),oma=c(5.1,4.1,4.1,2.1))
  }

#    x.lab <-  "Q" # not used I think
   xlim  <- range(pretty(BB))
   ylim  <- range(pretty(AA))
   #First the horizontal plot
   if(margin.plots) { 
      densF   <-  density(BB)
      plot(densF,ann=FALSE,xaxt="n",yaxt="n",type="l",xlim=xlim)   
      if(show.grid) grid()      
      title(ylab="Prob. Density",xpd=NA,mgp=c(1,1,0))
      title(main=maint,outer=TRUE)
      #Calculate 95% confidence intervals
      cumsumF.fun <-  approxfun(cumsum(densF$y)/sum(densF$y),densF$x)
      densF.fun   <-  approxfun(densF$x,densF$y)
      ul.F    <-  cumsumF.fun(1-alpha/2)      
      ul.dens <-  densF.fun(ul.F)
      ll.F    <-  cumsumF.fun(alpha/2)      
      ll.dens <-  densF.fun(ll.F)
      points(c(ll.F,ul.F),c(ll.dens,ul.dens),pch="|",cex=1.5)
      text(c(ll.F,ul.F),c(ll.dens,ul.dens),label=sprintf("%.3f",round(c(ll.F,ul.F),3)),pos=4)
      if(show.estimate) { 
        points(BB.est,densF.fun(BB.est),pch=19,cex=1.5)
        text(BB.est,densF.fun(BB.est),label=sprintf("%.3f",round(BB.est,3)),pos=4)
        }
      if(show.ll) { 
        points(LLb,densF.fun(LLb),pch=19,cex=1.5,col=4)
        text(LLb,densF.fun(LLb),label=sprintf("%.3f",round(LLb,3)),pos=4)
      }
    }
    #Now the vertical plot
    if(margin.plots) { 
      densAA <-  density(AA)
      plot(densAA$y,densAA$x,xaxt="n",yaxt="n",type="l",ylim=ylim)
      abline(v=0,col="grey")
      if(show.grid) grid()
      title(xlab="Prob. Density",xpd=NA,mgp=c(1,1,0))      
      #Calculate 95% confidence intervals
      cumsumAA.fun <-  approxfun(cumsum(densAA$y)/sum(densAA$y),densAA$x)
      densAA.fun   <-  approxfun(densAA$x,densAA$y)
      ul.AA    <-  cumsumAA.fun(1-alpha/2)      
      ul.dens <-  densAA.fun(ul.AA)
      ll.AA    <-  cumsumAA.fun(alpha/2)      
      ll.dens <-  densAA.fun(ll.AA)
      points(c(ll.dens,ul.dens),c(ll.AA,ul.AA),pch="-",cex=2)
      text(c(ll.dens,ul.dens),c(ll.AA,ul.AA),label=round(c(ll.AA,ul.AA),3),pos=4)
      if(show.estimate) {
        points(densAA.fun(AA.est),AA.est,pch=19,cex=1.5)      
        text(densAA.fun(AA.est),AA.est,,label=round(AA.est,3),pos=2)
        }
      if(show.ll) { 
        points(densAA.fun(LLa),LLa,pch=19,cex=1.5,col=4)      
        text(densAA.fun(LLa),LLa,,label=round(LLa,3),pos=2)
      }
    }
    #Now the main plot
    plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="")
#   Three sets of points for 3 Bayes MCMC chains - normally one would do
    if(show.points) points(BB[TH],AA[TH],pch=19,col=2,cex=0.05)
#    if(show.points) points(QC2[thin],QM2[thin]*0.15,pch=19,col=3,cex=0.05)
#    if(show.points) points(QC3[thin],QM3[thin]*0.15,pch=19,col=4,cex=0.05)
    title(xlab="B",ylab="A",xpd=NA)
#   Crude way to get second graph
    #if(show.points) points(QC1[thin],MZ1[thin],pch=19,col=2,cex=0.05)
    #if(show.points) points(QC2[thin],MZ2[thin],pch=19,col=3,cex=0.05)
    #if(show.points) points(QC3[thin],MZ3[thin],pch=19,col=4,cex=0.05)
    #title(xlab="Q",ylab="Z",xpd=NA)
    if(show.estimate) points(BB.est,AA.est,pch=19,cex=1.5)
    if(show.ll) points(LLb,LLa,pch=19,col=4,cex=1.5)
    if(show.grid) grid()
    contour(otolith.obj,levels=contour.args$levels,lables=NULL,add=TRUE)
#    lines below were origonal ---- I uses line above to get correct labled contours
# dont run lines below
    if(do.contours) {
      do.call(contour,otolith.obj)
    }
  #savePlot(filename = paste(my.dir,maint),type = "wmf") 

}

