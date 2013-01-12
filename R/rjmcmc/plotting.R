
############################################
# plotting routines
##############################################

LLplot <- function(fit) 
{
  lliks <- sapply(fit[1:3], function(x) sort(x $ llik))

  plot(0, 0, type = "n", 
       main = paste("LL of Bayes model:", fit $ stknam), 
       xlab = "model order", ylab = "log likelihood", 
       ylim = range(lliks), xlim = c(1, nrow(lliks)))
  for (i in 1:ncol(lliks)) 
  {
    lines(lliks[,i], lty = i, col = i)
  }
  lines(sort(fit $ SR $ llik), lwd = 2)
  legend(x = "bottomright", legend = c("HS","RK","BH","Avg"), 
         lty = c(1:ncol(lliks),1), col = c(1:ncol(lliks),1), lwd = c(1,1,1,2)) 
}


LLmodplot <- function(fit)
{
  rec <- fit $ data $ rec[order(fit $ data $ ssb)]
  ssb <- sort(fit $ data $ ssb)

  for (mod in 1:3) 
  { 
    model <- c("Hockey Stick", "Ricker", "Beverton Holt")[mod]
    Rsym <- exp(fit $ fits[[mod]])

    plot(ssb, rec, xlab = 'SSB', ylab = 'Recruits', 
         main = paste(fit $ stknam, model), 
         ylim = c(0, max(rec)), xlim = c(0, max(ssb)), type = "n")
    for (i in sample(1:nrow(Rsym), 1000)) {
      lines(ssb, Rsym[i,], col = paste0(grey(0), "10") )
    }
    for (i in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
      lines(ssb, apply(Rsym, 2, quantile, i), col = 2)
    }
    #TODO lines(ssb, Rsym[which.max(res[[mod]]),], col=1, lwd=2)
    lines(fit $ data $ ssb, fit $ data $ rec, col = "darkblue")
    points(ssb, rec, pch = 19, col = "darkblue")  
  }
}


# plot simulated recruits including Blim estimates
SRplot <- function (fit) 
{
  modset <- fit $ SR
  data <- fit $ data

  ssb <- data $ ssb
  rec <- data $ rec

  #TODO Blim <- res $ Blim
  #TODO PBlim <- res $ PBlim

  mn <- length(ssb)
  minSSB <- min(ssb, max(ssb)*0.05)
  maxSSB <- max(ssb)*1.1
  maxrec <- max(rec*1.5)

  plot(ssb, rec, xlim = c(0, maxSSB), ylim = c(0, maxrec), type = "n", 
       xlab = "SSB ('000 t)", ylab="Recruits", main = fit $ stknam)

  out <-
  do.call(rbind, lapply(1:1000, 
    function(i)
    {
      fssb <- runif(1000, minSSB, maxSSB)
      FUN <-  get(paste0(modset $ mod[i], "r"))
      frec <- exp( FUN(modset[i,], fssb) + rnorm(1000, sd = modset $ sigma[i]) )
      points(fssb, frec, pch = 20, col = paste0(grey(0), "10"), cex = 0.0625)

      data.frame(ssb = fssb, rec = frec)
    }))
  out $ grp <- with(out, floor(10 * (ssb - min(ssb)) / (max(ssb) - min(ssb) + 0.001)))
  out $ mid.grp <- with(out, (grp + 0.5) / 10 * (max(ssb) - min(ssb)) + min(ssb))

  #TODO use 
  summ <- with(out, 
    t(simplify2array( tapply(rec, grp, quantile, c(0.5, .025, .975)) )))

  mid.grp <- sort(unique(out $ mid.grp))

  lines(mid.grp, summ[,1], col = 7, lwd = 3)
  lines(mid.grp, summ[,2], col = 4, lwd = 3)
  lines(mid.grp, summ[,3], col = 4, lwd = 3)

  lines(ssb, rec, col = 10)
  points(ssb, rec, pch = 19, col = 10, cex = 1.25)


  #TODO plot Blim
  #lines(PBlim$mids,0.1*maxrec/max(PBlim$counts)*PBlim$counts,col=5,lwd=2)
  #lines(c(Blim,Blim),c(0,maxrec),col=5,lwd=2)
}




######## Creates equilibrium plots for a given Blim

Eqplot <- function (sim, stk, fit, Blim, Bpa = 1.4 * Blim)
{

  Nmod <- dim(sim $ ssbsa)[3]
  Nyrs <- dim(sim $ ssbsa)[2]

  Fscan <- sim $ Fscan

  catm <- apply(sim $ catsa, 1, mean)
  maxcatm <- which.max(catm)
  catsam <- apply(sim $ catsa, c(1,3), mean)
  maxpf <- apply(catsam, 2, which.max)
  fmsy <- Fscan[maxpf]
  v <- hist(fmsy, breaks = c(Fscan-Fscan[1]/2, max(Fscan) + Fscan[1]/2), plot = FALSE)

  pssb1 <- apply(sim $ ssbsa > Blim, 1, mean)
  pssb2 <- apply(sim $ ssbsa > Bpa, 1, mean)

  pp1 <- max(which(pssb1>.95))
  grad <- diff(Fscan[pp1 + 0:1]) / diff(pssb1[pp1 + 0:1])
  flim <- Fscan[pp1] + grad * (0.95 - pssb1[pp1])

  if (any(cumsum(v$counts) > Nmod / 2 - 1)) 
  {
    upp <- min(cumsum(v$counts)[ cumsum(v$counts) > Nmod / 2 - 1])
    uppp <- which(cumsum(v$counts)==upp)
  } else upp <- NA


  if (any(cumsum(v$counts) < Nmod / 2 - 1))
  {
    dwn <- max(cumsum(v$counts)[ cumsum(v$counts) < Nmod / 2 - 1])
    dwnp <- which(cumsum(v$counts)==dwn)
  } else dwn <- NA


  if (!is.na(upp) & !is.na(dwn)) 
  {
    vcum <- (Nyrs/2-dwn)*(Fscan[uppp]-Fscan[dwnp])/(upp-dwn)+Fscan[dwnp]   # 50%  on dist msy
  } else vcum <- NA

  vmode <- which.max(v$counts)                                                    # mode of dist msy
  msym <- sum(Fscan * v$counts) / sum(v$counts)                                   # mean of dist msy


  maint <- name(stk)  

  rec <- fit $ data $ rec
  ssb <- fit $ data $ ssb

  Catchs <- catch(stk)[, 1:length(ssb), drop = TRUE]
  FbarO <- fbar(stk)[, 1:length(ssb), drop = TRUE]

  recs <- sim $ recs
  ssbs <- sim $ ssbs
  cats <- sim $ cats
  ssbsa <- sim $ ssbsa

  NF <- length(Fscan)
  pp1 <- max(which(pssb1>.50))
  grad <- (Fscan[pp1+1]-Fscan[pp1])/(pssb3[pp1+1]-pssb3[pp1])
  flim5 <- Fscan[pp1]+grad*(0.5-pssb3[pp1])

  maxcatm <- which.max(catm)

  op <- par(mfrow = c(2, 2), mar = c(2.5, 2, 1.5, 1), oma = c(0, 0, 0, 0), 
            cex.axis = 0.75, tcl = 0.25, mgp = c(0, 0.25, 0), las = 1)

# recruits versus Fbar
  xmax <- max(Fscan)
  ymax <- max(recs[7,], rec)

  plot(Fscan, recs[7,], type = "l", lty = 4, 
       ylim = c(0, ymax), xlim = c(0, xmax), ylab = "", xlab = "")
  title(ylab = "Recruitment", xlab = "F bar", 
        cex.lab = 0.75, line = 1.2, cex.main = 0.75)
  mtext(text = paste(maint," a) Recruits"), cex = 0.75, side = 3, line = 0.5)
  lines(Fscan, recs[6,], lty = 3)
  lines(Fscan, recs[5,], lty = 2)
  lines(Fscan, recs[4,], lty = 1)
  lines(Fscan, recs[3,], lty = 2)
  lines(Fscan, recs[2,], lty = 3)
  lines(Fscan, recs[1,], lty = 4)
  points(FbarO, rec, pch = 21, cex = .75, bg = 1)
  lines(c(flim, flim), c(0, ymax), col = 3)

# recruits versus SSB
  ymax <- max(ssbs[7,])
  plot(Fscan, ssbs[7,], type = "l", lty = 4, 
       ylim = c(0, ymax), xlim = c(0, xmax), ylab = "", xlab = "")
  title(ylab = "SSB", xlab = "F bar", 
        cex.lab = 0.75, line = 1.2, cex.main = 0.75)
  mtext(text = "b) Spawning Stock Biomass", cex = 0.75, side = 3, line = 0.5)
  lines(Fscan, ssbs[6,], lty=3)
  lines(Fscan, ssbs[5,], lty=2)
  lines(Fscan, ssbs[4,], lty=1)
  lines(Fscan, ssbs[3,], lty=2)
  lines(Fscan, ssbs[2,], lty=3)
  lines(Fscan, ssbs[1,], lty=4)
  lines(c(0,xmax), c(Blim, Blim))
  text(x = 0.1, y = Blim * 1.1, "Blim", cex = 0.7)
  points(FbarO, ssb, pch = 21, cex = .75, bg = 1)
  lines(c(flim, flim), c(0, ymax), col = 3)


# catch versus Fbar
  ymax <- max(cats[7,])
  plot(Fscan, cats[7,], type = "l", lty = 4,
       ylim = c(0, ymax), xlim = c(0, max(Fscan)), ylab = "", xlab = "")
  title(ylab = "Catch", xlab = "F bar", 
        cex.lab = 0.75, line = 1.2, cex.main = 0.75)
  mtext(text = "c) Catch", cex = 0.75, side = 3, line = 0.5)
  #points(fbarsa[ssbsa<BLIM],catsa[ssbsa<BLIM],pch=21,col=6,bg=6,cex=0.125)
  #points(fbarsa[ssbsa>BLIM]+.0075,catsa[ssbsa>BLIM],pch=21,col=7,bg=7,cex=0.125)
  lines(Fscan, cats[7,], lty=4)
  lines(Fscan, cats[6,], lty=3)
  lines(Fscan, cats[5,], lty=2)
  lines(Fscan, cats[4,], lty=1)
  lines(Fscan, cats[3,], lty=2)
  lines(Fscan, cats[2,], lty=3)
  lines(Fscan, cats[1,], lty=4)
  points(FbarO, Catchs, pch = 21, cex = .75, bg = 1)
  lines(c(flim, flim), c(0, ymax), col = 3)
  lines(Fscan, catm, lty=1, col = 2)
  lines(rep(Fscan[maxcatm], 2), c(0, ymax), lty = 1, col = 5)

# F versus SSB
  plot(Fscan, 1-pssb1, type = "l", lty = 2,
       ylim = c(0,1), xlim = c(0,max(Fscan)), ylab = "", xlab = "")
  title(ylab = "Prob MSY, SSB<Bpa or Blim", xlab = "F bar", 
        cex.lab = 0.75, line = 1.2, cex.main = 0.75)
  mtext(text = "d) Prob MSY and Risk to SSB", cex = 0.75, side = 3, line = 0.5)
  lines(Fscan, 1 - pssb2, lty = 4)
  lines(c(0,xmax), c(0.05,0.05))
  lines(c(0,xmax), c(0.1,0.1))
  text(x = 0.1, y = 0.075, "5%", cex = 0.75)
  text(x = max(Fscan[pssb2 > 0.5]) - .05, y = 0.5, "SSB<Bpa", cex = 0.75)
  text(x = max(Fscan[pssb1 > 0.7]) + .1, y = 0.3, "SSB<Blim", cex = 0.75)
  lines(c(flim,flim), c(0,1), col = 3)
  lines(v$mids, v$counts / Nmod, col = 4)
  text(x = 0.2, y = 0.15, "Prob of Fmsy", cex = 0.75)
  lines(rep(vcum,2), c(0,1), lty = 1, col = 5)

  FCrash5  <- Fscan[ which(cats[2, which.max(cats[2,]):NF] < 0.5*max(cats[2,]) )[1] ]
  FCrash50 <- Fscan[ which(cats[4, which.max(cats[4,]):NF] < 0.5*max(cats[4,]) )[1] ]

  out <- c(Blim, Bpa, flim, flim5, round(vcum*100)/100, Fscan[maxcatm], FCrash5, FCrash50)
  names(out) <- c("Blim","Bpa","Flim","Flim5","MSY:median","Maxmeanland","FCrash5","FCrash50")

  out
}

