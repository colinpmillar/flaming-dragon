

library(FLCore)
library(colorout)
library(numDeriv)
library(nlstools)
library(scam)
library(xtable)

source('R/functions.R')
  
data(ple4)
stk <- ple4
ices.Fmsy <- 0.3 # no always there for our stocks

### Ref point biological data ###;

yr <- "2008"

matp         <- mat(stk)[, yr, drop=TRUE] 
mp           <- m(stk)[, yr, drop=TRUE]
age          <- as.numeric(names(matp))
n.age        <- length(age)
stock.weight <- stock.wt(stk)[, yr, drop=TRUE]
catch.weight <- catch.wt(stk)[, yr, drop=TRUE]
Sel          <- harvest(stk)[, yr, drop=TRUE]
pr           <- Sel/fbar(stk)[, yr, drop=TRUE]


########################################################;

stock.name <- name(stk) 
rec.age <- min(age) # might want to change this
sscale <- 1000
s.name <- c('SSB (Kt)') # check units
rscale <- 1000
r.name <- c('  Recruit (millions)') # check units
p.name <- c('Recruit per spawner  ')


n <- with(dims(stk), maxyear - minyear + 1)
s <- ssb(stk)[, 1:(n-rec.age), drop=TRUE]
r <- stock.n(stk)[paste(rec.age), (rec.age+1):n, drop=TRUE]
year <- as.numeric(names(r))

filt <- !is.na(r) | !is.na(s)
r    <- r[filt]
s    <- s[filt]
year <- year[filt]

start.alpha.bh <- max(r)
ms <- mean(s)
mr <- mean(r)
n  <- length(r)

x <- s/sscale
y <- r/rscale
ys <- y/x
log.ys <- log(ys)   

wt.ext <- rep(1, n)  

srdat <- data.frame(logrec = log(y), logssb = log(x), ssb = x, weights = wt.ext)

### Beverton-Holt  #####

init.Rmax <- 1000
init.S50 <- 300
init.alpha <- init.Rmax
init.beta <- init.S50
ls.slope <- sum(x * y) / sum(x * x) 
Rmax <- max(y)
Pmax <- 2 * max(y / x)

BH.fit <- nls(logrec ~ log(alpha) + logssb - log(beta + ssb), data = srdat,
           start = list(alpha = init.alpha, beta = init.beta),
           lower = c(0.1,1), upper = c(Rmax,1000),
           algorithm = "port")

summary(BH.fit)
BH.parm <- coef(BH.fit)

ssb.max.bh <- max(x) * 6 # to give some room
ssb.pred.points.bh <- seq(0.1, ssb.max.bh, length = 500)
BH.pred.rec <- exp( predict(BH.fit, list(ssb = ssb.pred.points.bh, logssb = log(ssb.pred.points.bh))) )  
BH.pred.prod <- BH.pred.rec / ssb.pred.points.bh

BH.raw.res <- resid(BH.fit)
BH.sig2 <- sum(BH.raw.res^2) / (sum(wt.ext) - 2)
BH.res <- BH.raw.res / sqrt(BH.sig2)

BH.boo <- nlsBoot(BH.fit)

pdf(file = 'fig/BHparms_boo.pdf')
plot(BH.boo)
dev.off()

ssb.max <- max(x) * 3
ssb.pred.points <- seq(0.1, ssb.max, length = 500)

##### Get BH MSY stuff  #######;

bh.ref.points <- make.RP(log(BH.parm), 0.2, sr_type = "BH")
BH.Bmsy <- bh.ref.points[4]

BH.boo.RP <- apply(log(BH.boo $ coef), 1, make.RP, Fstart = 0.1, sr_type = "BH")
BH.boo.qRP <- apply(BH.boo.RP,1,quantile,probs=c(0.025,0.25,0.5,0.75,0.975))

### Hockey-Stick ##########

init.alpha = 0.5 * init.Rmax / init.S50
init.delta = init.S50
gam2by4 = 0.1 / 4
delta.max = max(x) 
delta.min = min(x)

#rec = log(alpha) + log(ssb + sqrt(delta**2 + gam2by4) - sqrt((ssb-delta)**2 + gam2by4))
#plot(ssb,rec)

HS.fit <- nls(logrec ~ log(alpha) + log(ssb + sqrt(delta^2 + gam2by4) - sqrt((ssb-delta)^2 + gam2by4)),
           data=srdat,
           start = list(alpha = init.alpha,delta = init.delta),
           lower = c(0, delta.min), upper = c(25, delta.max),
           algorithm = "port")

summary(HS.fit)
HS.parm = coef(HS.fit)

HS.pred.rec = exp(predict(HS.fit,list(ssb=ssb.pred.points)))  
HS.pred.prod =  HS.pred.rec/ssb.pred.points

HS.raw.res = resid(HS.fit)
HS.sig2 = sum(HS.raw.res^2) / (sum(wt.ext)-2)
HS.res = HS.raw.res/sqrt(HS.sig2)

HS.boo <- nlsBoot(HS.fit)

pdf(file = 'fig/HSparms_boo.pdf')
plot(HS.boo)
dev.off()

##### Get HS MSY stuff  #######;

hs.ref.points = make.RP(log(HS.parm), 0, sr_type = "HS")
HS.Bmsy = hs.ref.points[4]

HS.boo.RP = apply(log(HS.boo$coefboot),1,make.RP,Fstart=0, sr_type = "HS");
HS.boo.qRP = apply(HS.boo.RP,1,quantile,probs=c(0.025,0.25,0.5,0.75,0.975))

### Ricker ##########

init.alpha= init.Rmax
init.lbeta = -5
srdat$one = 1

RK.fit <- nls(logrec ~ log(alpha) + lbeta + one + logssb - exp(lbeta)*ssb,
           data = srdat,
           start = list(alpha = init.alpha, lbeta = init.lbeta),
           lower = c(0, -Inf), upper = c(Rmax, 0),
           algorithm = "port")

summary(RK.fit)
RK.parm = coef(RK.fit)

RK.pred.rec = exp(predict(RK.fit,list(ssb=ssb.pred.points,logssb=log(ssb.pred.points), one = rep(1,length(ssb.pred.points)))))  
RK.pred.prod =  RK.pred.rec/ssb.pred.points

#RK.mean.rec = exp(predict(RK.fit,list(ssb=mx,logssb=log(mx),one = 1)))

RK.raw.res = resid(RK.fit)
RK.sig2 = sum(RK.raw.res^2)/(sum(wt.ext)-2)
RK.res = RK.raw.res/sqrt(RK.sig2)

RK.boo <- nlsBoot(RK.fit)

pdf(file = 'fig/RKparms_boo.pdf')
plot(RK.boo)
dev.off()

##### Get RK MSY stuff  #######;


RK.parmi = c(RK.parm[1]*(exp(1)*exp(RK.parm[2])),exp(RK.parm[2]))
rk.ref.points = make.RP(log(RK.parmi), 0.2, sr_type = 'RK')
RK.Bmsy = rk.ref.points[4]

RK.boo$coefbooti = RK.boo$coefboot
RK.boo$coefbooti[,2] = exp(RK.boo$coefbooti[,2])
RK.boo$coefbooti[,1] = RK.boo$coefboot[,1] * (exp(1) * RK.boo$coefbooti[,2])

RK.boo.RP = apply(log(RK.boo$coefbooti),1,make.RP, Fstart=0.1, sr_type = 'RK');
RK.boo.qRP = apply(RK.boo.RP,1,quantile,probs=c(0.025,0.25,0.5,0.75,0.975))

####  shape constrained Spline fits #########
    
rng = diff(range(x)) / n
ssb.pred.points1 = c(seq(0.1, min(x), by = rng/4), seq(min(x), ssb.max, by=rng))
    
np = length(ssb.pred.points1)
dat = data.frame(stock.size = c(x,ssb.pred.points1), recruit = c(y,rep(1,np)), wt=c(wt.ext,rep(0,np)))
dat $ log.recruit = log(dat$recruit)
dat $ offset = log(dat$stock.size)

## Need this code to constrain prod at low stock size;
dat1 = dat
pos = which.min(dat1$stock.size)
dat1 $ recruit[pos] = Pmax*dat1$stock.size[pos]
dat1 $ log.recruit[pos] = log(dat1$recruit[pos])  
dat1 $ wt[pos] = 10000

dato = dat

dat.nz = subset(dat,wt==1)

## Nonparametric compensatory mortality SR model;

m = length(dat$stock.size)
nknots = 20

mygcv <- function(rho) 
{
  temp <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd",m=2) + offset(offset),
          family=gaussian(link="identity"),data=dat,optimizer="nlm",
          weights=dat$wt,sp=exp(rho))
  c(n * temp $ dev / (n - sum(temp $ edf))^2 )
}

gcv.fit = optim(-7, mygcv, upper=5, lower=-8, method="L-BFGS-B")

b.cm <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd",m=2) + offset(offset),
        family=gaussian(link="identity"),data=dat,optimizer="nlm",
        weights=dat$wt,sp=exp(gcv.fit$par))
       
ind = order(dat$stock.size)
plot(dat$stock.size[ind],exp(b.cm$fitted.values[ind]),type='l',ylim=c(0,max(y)))
points(x,y)
prod = exp(b.cm$fitted.values[ind])/dat$stock.size[ind]
plot(dat$stock.size[ind],prod,type='l',ylim=c(0,Pmax))
points(x,y/x)

#prod1 = exp(b.cm1$fitted.values[ind])/dat$stock.size[ind]
#lines(dat$stock.size[ind],prod1,col='red')
                
#dat=dato

bootscam = function(scamfit,sp,nboot){
  ret=matrix(NA,length(dat$stock.size),nboot)
  ind = dat$wt==1
  raw.resid = scamfit$residuals[ind]
  for(i in 1:nboot){ 
    boot.dat=dat
    boot.y = exp(log(y) + sample(raw.resid,n,replace=T))
    boot.dat$recruit = c(boot.y,rep(1,np))
    boot.dat$log.recruit=log(boot.dat$recruit)
    scamfiti = scam(scamfit$formula,family=scamfit$family,optimizer=scamfit$optimizer,
               dat=boot.dat,weights=boot.dat$wt,sp=sp,start=scamfit$coefficients)
    pred.rec = exp(scamfiti$fitted.values)      
    pred.prod = pred.rec/boot.dat$stock.size
    max.pred.rec = max(pred.rec)      
    max.pred.prod = max(pred.prod)

#    if(max.pred.rec>Rmax){
#      npp = length(boot.dat$recruit)
#      boot.dat$recruit[npp] = Rmax 
#      boot.dat$log.recruit[npp]=log(Rmax)  
#      boot.dat$wt[npp]=10000    
#      scamfiti = scam(scamfit$formula,family=scamfit$family,optimizer=scamfit$optimizer,
#                 dat=boot.dat,weights=boot.dat$wt,sp=sp)
#     }
#     if(max.pred.prod>Pmax){
#      pos = which.min(boot.dat$stock.size)
#      boot.dat$recruit[pos] = Pmax*boot.dat$stock.size[pos]
#      boot.dat$log.recruit[pos]=log(boot.dat$recruit[pos])  
#      boot.dat$wt[pos]=10000
#      scamfiti = scam(scamfit$formula,family=scamfit$family,optimizer=scamfit$optimizer,
#                 dat=boot.dat,weights=boot.dat$wt,sp=sp)
#     }   
    ret[,i] = scamfiti$fitted.values 
  }      
return(ret)
}

nboot=2000
boot.b.cm = exp(bootscam(b.cm,exp(gcv.fit$par),nboot))
boot.b.cm = exp(bootscam(b.cm,exp(gcv.fit$par),nboot))
boot.b.cmq = apply(boot.b.cm,1,quantile,probs=c(0.025,0.975))
       
mydf <- function(rho){
temp <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd",m=2) + offset(offset),
        family=gaussian(link="identity"),data=dat,optimizer="nlm",
        weights=dat$wt,sp=exp(rho))
ret = (sum(temp$edf) - my.df)**2                                  
return(ret)
}

my.df = 3
fit3 = optim(-7,mydf,upper=5,method="L-BFGS-B")
sp3 = exp(fit3$par)
b.cm3 <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd",m=2) + offset(offset),
        family=gaussian(link="identity"),data=dat,optimizer="nlm",
        weights=dat$wt,sp=sp3)
        
boot.b.cm3 = exp(bootscam(b.cm3,sp3,nboot))
boot.b.cm3q = apply(boot.b.cm3,1,quantile,probs=c(0.025,0.975))
        
my.df = 4
fit4 = optim(-7,mydf,upper=5,method="L-BFGS-B")
sp4 = exp(fit4$par)
b.cm4 <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd",m=2) + offset(offset),
        family=gaussian(link="identity"),data=dat,optimizer="nlm",
        weights=dat$wt,sp=sp4)

boot.b.cm4 = exp(bootscam(b.cm4,sp4,nboot))
boot.b.cm4q = apply(boot.b.cm4,1,quantile,probs=c(0.025,0.975))

## scam does not adjust for zero-weighted data when computing df's;
se.scale = sqrt((np+n-sum(b.cm$edf))/(n-sum(b.cm$edf)))
b.cm$fitted.se = se.scale*sqrt(diag(b.cm$X %*% b.cm$Vp.t %*% t(b.cm$X)))
b.cm$pred.prod = exp(b.cm$fitted.value)/dat$stock.size

se.scale = sqrt((np+n-sum(b.cm3$edf))/(n-sum(b.cm3$edf)))
b.cm3$fitted.se = se.scale*sqrt(diag(b.cm3$X %*% b.cm3$Vp.t %*% t(b.cm3$X)))
b.cm3$pred.prod = exp(b.cm3$fitted.value)/dat$stock.size

se.scale = sqrt((np+n-sum(b.cm4$edf))/(n-sum(b.cm4$edf)))
b.cm4$fitted.se = se.scale*sqrt(diag(b.cm4$X %*% b.cm4$Vp.t %*% t(b.cm4$X)))
b.cm4$pred.prod = exp(b.cm4$fitted.value)/dat$stock.size

            
#summary(b.ccm)


##### Stock-Recruit Plot ##### 

ind = order(dat$stock.size)
  
b.cm.pred = data.frame(stock.size = dat$stock.size[ind],recruit = exp(b.cm$fitted.values)[ind],
  LCI.rec.asy = exp(b.cm$fitted.values - qnorm(0.975)*b.cm$fitted.se)[ind],
  UCI.rec.asy = exp(b.cm$fitted.values + qnorm(0.975)*b.cm$fitted.se)[ind],    
  LCI.rec = boot.b.cmq[1,ind],
  UCI.rec = boot.b.cmq[2,ind],
  prod = b.cm$pred.prod[ind],
  recruit4 = exp(b.cm4$fitted.values)[ind],     
  LCI.rec4 = boot.b.cm4q[1,ind],
  UCI.rec4 = boot.b.cm4q[2,ind])
b.cm.pred$LCI.prod = b.cm.pred$LCI.rec/b.cm.pred$stock.size  
b.cm.pred$UCI.prod = b.cm.pred$UCI.rec/b.cm.pred$stock.size  
b.cm.pred$prod4 = b.cm.pred$recruit4/b.cm.pred$stock.size   
b.cm.pred$LCI.prod4 = b.cm.pred$LCI.rec4/b.cm.pred$stock.size  
b.cm.pred$UCI.prod4 = b.cm.pred$UCI.rec4/b.cm.pred$stock.size 
  
b.cm.predo = subset(b.cm.pred,((stock.size<=max(x))&(stock.size>=min(x))))
  

for (c in 1:2){

if(c==1){xlim = range(b.cm.pred$stock.size)}  
if(c==2){xlim = range(b.cm.predo$stock.size)}

file.name <- paste(path.to.file,'fit',c,'.png', sep = "") 
png(file=file.name,width=3,height=3,units='in',res=1200) 
par(mfrow=c(2,1),oma=c(3,1,0.5,0),mar=c(0.5,3,0,1),las=0,cex=0.85)

ylim1 = c(0,max(b.cm.pred$UCI.rec,b.cm.pred$UCI.rec4,BH.pred.rec,RK.pred.rec,y))
ylim2 = c(0,max(b.cm.predo$UCI.rec,b.cm.predo$UCI.rec4,BH.pred.rec,RK.pred.rec,y))
if(ylim1[2] > 1.5*ylim2[2]){ylim1[2]=1.5*ylim2[2]}
ylim1o=ylim1
if(c==1){ylim = ylim1}  
if(c==2){ylim = ylim2}  
     
plot(b.cm.pred$stock.size,b.cm.pred$recruit,type='l',lwd=1.5,
  ylab='',xlab='',las=1,xaxt='n',ylim=ylim,xlim=xlim)  
points(x,y,cex=0.5)  
lines(b.cm.pred$stock.size,b.cm.pred$LCI.rec,type='l',lwd=1,lty=2)
lines(b.cm.pred$stock.size,b.cm.pred$UCI.rec,type='l',lwd=1,lty=2)
#lines(ssb.pred.points,HS.pred.rec,lty=1,col='green',lwd=1.5)
lines(ssb.pred.points.bh,BH.pred.rec,lty=1,col='blue',lwd=1.5)  
lines(ssb.pred.points,RK.pred.rec,lty=1,col='red',lwd=1.5)   
lines(b.cm.pred$stock.size,b.cm.pred$recruit4,lty=1,col='purple',lwd=1.5) 
lines(b.cm.pred$stock.size,b.cm.pred$LCI.rec4,type='l',lwd=1,lty=2,col='purple')
lines(b.cm.pred$stock.size,b.cm.pred$UCI.rec4,type='l',lwd=1,lty=2,col='purple')  
 
mtext(side=2,line=3,outer=FALSE,r.name,cex=0.85)

ind1 = ((ssb.pred.points<=max(x))&(ssb.pred.points>=min(x)))
ylim1 = c(0,max(b.cm.pred$UCI.prod,HS.pred.prod,BH.pred.prod,RK.pred.prod,y/x))
ylim2 = c(0,max(b.cm.predo$UCI.prod,HS.pred.prod[ind1],BH.pred.prod[ind1],RK.pred.prod[ind1],y/x))
if(ylim1[2] > 1.5*ylim2[2]){ylim1[2]=1.5*ylim2[2]}
if(c==1){ylim = ylim1}  
if(c==2){ylim = ylim2}
       
plot(b.cm.pred$stock.size,b.cm.pred$prod,type='l',lwd=1.5,
  ylab='',xlab='',las=1,ylim=ylim,xlim=xlim)        
points(x,y/x,cex=0.5)    
lines(b.cm.pred$stock.size,b.cm.pred$LCI.prod,type='l',lwd=1,lty=2)
lines(b.cm.pred$stock.size,b.cm.pred$UCI.prod,type='l',lwd=1,lty=2)
#lines(ssb.pred.points,HS.pred.prod,lty=1,col='green',lwd=1.5)
lines(ssb.pred.points.bh,BH.pred.prod,lty=1,col='blue',lwd=1.5) 
lines(ssb.pred.points,RK.pred.prod,lty=1,col='red',lwd=1.5)  
lines(b.cm.pred$stock.size,b.cm.pred$prod4,lty=1,col='purple',lwd=1.5) 
lines(b.cm.pred$stock.size,b.cm.pred$LCI.prod4,type='l',lwd=1,lty=2,col='purple')
lines(b.cm.pred$stock.size,b.cm.pred$UCI.prod4,type='l',lwd=1,lty=2,col='purple')      

legend("topright",lty=1,col=c('black','purple','blue','red'),lwd=1.5,
c('NP GCV','NP df=4','BH','RK'),bty='n',cex=0.85) 
  
mtext(side=2,line=3,outer=FALSE,p.name,cex=0.85)     
mtext(side=1,line=2.2,outer=FALSE,s.name,cex=0.85)

dev.off()
}

#### compare smoothers with different df's #####

ind = order(dat$stock.size)

b.cm.comp = data.frame(stock.size = dat$stock.size[ind],
  recruit = exp(b.cm$fitted.values)[ind],     
  LCI.rec = boot.b.cmq[1,ind],
  UCI.rec = boot.b.cmq[2,ind],
  recruit3 = exp(b.cm3$fitted.values)[ind],     
  LCI.rec3 = boot.b.cm3q[1,ind],
  UCI.rec3 = boot.b.cm3q[2,ind],
  recruit4 = exp(b.cm4$fitted.values)[ind],     
  LCI.rec4 = boot.b.cm4q[1,ind],
  UCI.rec4 = boot.b.cm4q[2,ind])
  
par(mar=c(4,4,0.5,1),las=0)

#ylim = c(0,max(b.cm.comp$UCI.rec,b.cm.comp$UCI.rec3,b.cm.comp$UCI.rec4))
#ylim[2]=10 
     
plot(b.cm.comp$stock.size,b.cm.comp$recruit,type='l',lwd=2,
  ylab='',xlab='',las=1,ylim=ylim1o)
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec,type='l',lwd=1,lty=2)
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec,type='l',lwd=1,lty=2)
    
lines(b.cm.comp$stock.size,b.cm.comp$recruit3,type='l',lwd=2,lty=1,col='brown')
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec3,type='l',lwd=1,lty=2,col='brown')
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec3,type='l',lwd=1,lty=2,col='brown') 
    
lines(b.cm.comp$stock.size,b.cm.comp$recruit4,type='l',lwd=2,lty=1,col='purple')
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec4,type='l',lwd=1,lty=2,col='purple')
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec4,type='l',lwd=1,lty=2,col='purple')

points(x,y)
    
legend("topright",lty=1,lwd=2,c('GCV','df=3','df=4'),bty='n',col=c('black','brown','purple'))
 
mtext(side=2,line=3,outer=FALSE,r.name)     
mtext(side=1,line=2.2,outer=FALSE,s.name)

file.name <- paste(path.to.file,'sp_compare', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)
dev.off()


#### compare smoothers with different df's #####

ind = order(dat$stock.size)

b.cm.comp = data.frame(stock.size = dat$stock.size[ind],
  recruit = exp(b.cm$fitted.values)[ind],     
  LCI.rec = boot.b.cmq[1,ind],
  UCI.rec = boot.b.cmq[2,ind],
  recruit3 = exp(b.cm3$fitted.values)[ind],     
  LCI.rec3 = boot.b.cm3q[1,ind],
  UCI.rec3 = boot.b.cm3q[2,ind],
  recruit4 = exp(b.cm4$fitted.values)[ind],     
  LCI.rec4 = boot.b.cm4q[1,ind],
  UCI.rec4 = boot.b.cm4q[2,ind])
  
par(mar=c(4,4,0.5,1),las=0)

#ylim = c(0,max(b.cm.comp$UCI.rec,b.cm.comp$UCI.rec3,b.cm.comp$UCI.rec4))
#ylim[2]=10 
     
plot(b.cm.comp$stock.size,b.cm.comp$recruit,type='l',lwd=2,
  ylab='',xlab='',las=1,ylim=ylim1o)
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec,type='l',lwd=1,lty=2)
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec,type='l',lwd=1,lty=2)
    
lines(b.cm.comp$stock.size,b.cm.comp$recruit3,type='l',lwd=2,lty=1,col='brown')
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec3,type='l',lwd=1,lty=2,col='brown')
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec3,type='l',lwd=1,lty=2,col='brown') 
    
lines(b.cm.comp$stock.size,b.cm.comp$recruit4,type='l',lwd=2,lty=1,col='purple')
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec4,type='l',lwd=1,lty=2,col='purple')
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec4,type='l',lwd=1,lty=2,col='purple')

points(x,y)
    
legend("right",lty=1,lwd=2,c('GCV','df=3','df=4'),bty='n',col=c('black','brown','purple'))
 
mtext(side=2,line=3,outer=FALSE,r.name)     
mtext(side=1,line=2.2,outer=FALSE,s.name)

file.name <- paste(path.to.file,'sp_compare', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)
dev.off()

##### Residuals Plots #####

resid.matrix = matrix(NA,n,4)
smooth.resid.matrix = resid.matrix
ind = dat$wt==1 
sort.x = sort(x)

resid.matrix[,1] = b.cm$residuals[ind] 
resid.matrix[,2] = HS.raw.res      
resid.matrix[,3] = BH.raw.res      
resid.matrix[,4] = RK.raw.res

resid.scale = sqrt(mean(resid.matrix**2)) 

smooth.resid.matrix[,1] = predict(loess(resid.matrix[,1]~x,span=0.5),sort.x)/resid.scale
smooth.resid.matrix[,2] = predict(loess(resid.matrix[,2]~x,span=0.5),sort.x)/resid.scale
smooth.resid.matrix[,3] = predict(loess(resid.matrix[,3]~x,span=0.5),sort.x)/resid.scale  
smooth.resid.matrix[,4] = predict(loess(resid.matrix[,4]~x,span=0.5),sort.x)/resid.scale

ylim = c(min(smooth.resid.matrix),max(smooth.resid.matrix))
    
par(mar=c(4,4,0.5,1),las=0)
       
plot(sort.x,smooth.resid.matrix[,1],type='l',ylab='',xlab='',las=1,col='black',ylim=ylim,lwd=2)
lines(sort.x,smooth.resid.matrix[,2],col='green',lwd=2)                                       
lines(sort.x,smooth.resid.matrix[,3],col='blue',lwd=2)                                        
lines(sort.x,smooth.resid.matrix[,4],col='red',lwd=2)

abline(h=0,lty=2,lwd=2)
abline(h=-0.5,lty=2,lwd=2,col='grey')
abline(h=0.5,lty=2,lwd=2,col='grey')

legend("topright",lty=1,col=c('black','green','blue','red'),lwd=2,
c('CM','HS','BH','RK'),bty='n')
       
mtext(side=2,line=3,"Standardized residual pattern")
mtext(side=1,line=2.5,s.name)

file.name <- paste(path.to.file,'resid', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)

dev.off()


#####  MSY Plots  ############ 

xsect <- function(diff,ssb){
  if(all(diff<0)){pssb=0}
  if(any(diff>0)){
    ipos = 1:length(diff)
    p1 = min(ipos[diff<0])
    p2 = max(ipos[diff>=0])
    wt = 1 - diff[p2]/(diff[p2]-diff[p1])
    pssb = wt*ssb[p2] + (1-wt)*ssb[p1]
  }
  return(pssb)
}

scam.equil = function(r,s,F){
  temp1 = c(1/BPR,1,YPR/BPR)
  diff = r - s/BPR
  ret = c(0,0,0)
  if(any(diff<0)){ret = temp1*xsect(diff,s)}    
  if(all(diff<0)){ret = temp1*NA}      
  if(all(diff>0)){ret = Rmax*BPR}
return(ret)
}

scam.rp = function(F){
  fp = F*pr                
  zp = mp + fp
  czp = cumsum(c(0,zp[1:(n.age-1)]))
  YPR = sum(exp(-czp)*catch.weight*fp*(1-exp(-zp))/zp)
  BPR = sum(exp(-czp - zp/2)*stock.weight*matp)
  temp1 = c(1/BPR,1,YPR/BPR)
  diff = r - s/BPR
  ret = c(0,0,0)
  if(any(diff<0)){ret = temp1*xsect(diff,s)} 
  if(all(diff<0)){ret = temp1*NA}      
  if(all(diff>0)){ret = Rmax*BPR*temp1}
  ret
}

scam.msy = function(F){-scam.rp(F)[3]}
   
Fv = seq(0.1,1.0,by=0.005)
nFv = length(Fv) 
mFv = matrix(Fv,nrow=nFv,1)

BH.equil = matrix(NA,nrow=nFv,ncol=3) 
RK.equil = matrix(NA,nrow=nFv,ncol=3) 
HS.equil = matrix(NA,nrow=nFv,ncol=3)
CM.equil = matrix(NA,nrow=nFv,ncol=3)
CM3.equil = matrix(NA,nrow=nFv,ncol=3)
CM4.equil = matrix(NA,nrow=nFv,ncol=3)

for(i in 1:nFv){
  F = Fv[i]
  sr_type='BH'
  temp = make.RP(log(BH.parm),F,0)
  if(temp[4]<ssb.max.bh){BH.equil[i,] = temp[3:5]} 
  sr_type='RK'
  temp = make.RP(log(RK.parmi),F,0)
  if(temp[4]<ssb.max){RK.equil[i,] = temp[3:5]} 
  sr_type='HS'
  temp = make.RP(log(HS.parm),F,0)
  if(temp[4]<ssb.max){HS.equil[i,] = temp[3:5]}

  fp = F*pr                
  zp = mp + fp
  czp = cumsum(c(0,zp[1:(n.age-1)]))
  YPR = sum(exp(-czp)*catch.weight*fp*(1-exp(-zp))/zp)
  BPR = sum(exp(-czp - zp/2)*stock.weight*matp)

#ylim = range(b.cm.pred$recruit)
#plot(ssb.pred.points,BH.pred.rec,lty=1,col='blue',lwd=2,type='l',ylim=ylim) 
#lines(ssb.pred.points,HS.pred.rec,lty=1,col='green',lwd=2) 
#lines(b.cm.pred$stock.size,b.cm.pred$recruit,lty=1,col='black',lwd=2)
#abline(a=0,b=1/BPR,lty=1)
  
   CM.equil[i,] = scam.equil(b.cm.comp$recruit,b.cm.comp$stock.size,F)
   CM3.equil[i,] = scam.equil(b.cm.comp$recruit3,b.cm.comp$stock.size,F)
   CM4.equil[i,] = scam.equil(b.cm.comp$recruit4,b.cm.comp$stock.size,F)
    
}

par(mar=c(5,5,0.5,0.5))
ylim = range(CM.equil[,3],CM3.equil[,3],CM4.equil[,3],na.rm=T)
plot(Fv,CM.equil[,3],ylim=ylim,lwd=2,col='black',type='l',xlab='Fishing mortality',ylab='Yield (Kt)',las=1)
lines(Fv,CM3.equil[,3],lwd=2,lty=1,col='brown')     
lines(Fv,CM4.equil[,3],lwd=2,lty=1,col='purple')     
legend("bottom",lty=1,col=c('black','brown','purple'),lwd=2,
c('GCV','df=3','df=4'),bty='n')

file.name <- paste(path.to.file,'yield1', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)

dev.off()

file.name <- paste(path.to.file,'yield.png', sep = "")
png(file=file.name,width=2.5,height=2.5,units='in',res=1200) 
par(mar=c(3.5,3.5,0.5,0.5),las=0,cex=0.85)

ylim = range(BH.equil[,3],RK.equil[,3],HS.equil[,3],CM.equil[,3],CM4.equil[,3],na.rm=T)
xlim=c(0.1,1.0)
plot(Fv,BH.equil[,3],ylim=ylim,lwd=2,col='blue',type='l',xlab='',ylab='',las=1,xlim=xlim)
lines(Fv,HS.equil[,3],lwd=2,lty=1,col='green')     
lines(Fv,CM.equil[,3],lwd=2,lty=1,col='black')       
lines(Fv,RK.equil[,3],lwd=2,lty=1,col='red')         
lines(Fv,CM4.equil[,3],lwd=2,lty=1,col='purple')       
legend("bottomright",lty=1,col=c('black','purple','green','blue','red'),lwd=2,
 c('NP GCV','NP df=4','HS','BH','RK'),bty='n',cex=0.8,inset=c(0.1,0)) 
mtext(side=2,line=2.5,"Fishing mortality")
mtext(side=1,line=2.5,"Yield (Kt)")
  
dev.off()

# get scam bootstrapped MSY RP's;

Fv = seq(0.1,3,by=0.005)
nFv = length(Fv) 
mFv = matrix(Fv,nrow=nFv,1)

CM.boo.RP = matrix(NA,3,nboot);
CM3.boo.RP = matrix(NA,3,nboot);
CM4.boo.RP = matrix(NA,3,nboot);

ind = order(dat$stock.size)
s = dat$stock.size[ind]
for(i in 1:nboot){

  r=boot.b.cm[ind,i]
  tf = apply(mFv,1,scam.msy)
  tf[tf==0]=NA  
  if(any(!is.na(tf))){
    start.f=Fv[!is.na(tf)]
    fs.pos = which.min(tf)
    fs=mFv[fs.pos];fl=min(start.f);fu=max(start.f)
    if(fu>fl){
      fit = optim(fs,scam.msy,upper=fu,lower=fl,method="L-BFGS-B")
      CM.boo.RP[,i]=scam.rp(fit$par);CM.boo.RP[1,i]=fit$par
    }
    if(fl==fu){CM.boo.RP[,i]=scam.rp(fu);CM.boo.RP[1,i]=fu}
  }
  
  r=boot.b.cm3[ind,i]
  tf = apply(mFv,1,scam.msy)
  tf[tf==0]=NA
  if(any(!is.na(tf))){
    start.f=Fv[!is.na(tf)]
    fs.pos = which.min(tf)
    fs=mFv[fs.pos];fl=min(start.f);fu=max(start.f)
    if(fu>fl){
      fit = optim(fs,scam.msy,upper=fu,lower=fl,method="L-BFGS-B")
      CM3.boo.RP[,i]=scam.rp(fit$par);CM3.boo.RP[1,i]=fit$par
    }
    if(fl==fu){CM3.boo.RP[,i]=scam.rp(fu);CM3.boo.RP[1,i]=fu}
  } 
  
  r=boot.b.cm4[ind,i]
  tf = apply(mFv,1,scam.msy)
  tf[tf==0]=NA
  if(any(!is.na(tf))){
    start.f=Fv[!is.na(tf)]
    fs.pos = which.min(tf)
    fs=mFv[fs.pos];fl=min(start.f);fu=max(start.f)
    if(fu>fl){
      fit = optim(fs,scam.msy,upper=fu,lower=fl,method="L-BFGS-B")
      CM4.boo.RP[,i]=scam.rp(fit$par);CM4.boo.RP[1,i]=fit$par
    }
    if(fl==fu){CM4.boo.RP[,i]=scam.rp(fu);CM4.boo.RP[1,i]=fu}
  }
 # print(i)
}

CM.boo.qRP = apply(CM.boo.RP,1,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
colnames(CM.boo.qRP) = c('Fmsy','Bmsy','MSY')
CM3.boo.qRP = apply(CM3.boo.RP,1,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
colnames(CM3.boo.qRP) = c('Fmsy','Bmsy','MSY')
CM4.boo.qRP = apply(CM4.boo.RP,1,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
colnames(CM4.boo.qRP) = c('Fmsy','Bmsy','MSY')


##### MSY table #####

YPR <- function(f){
    n.age = length(pr)
    fp = f*pr
    zp = mp + fp
    czp = cumsum(c(0,zp[1:(n.age-1)]))
    YPR = sum(exp(-czp)*catch.weight*fp*(1-exp(-zp))/zp)
    return(YPR)
}   

ypr3 = apply(mFv,1,YPR)
Fmax = mFv[which.max(ypr3)]

msy.tab = matrix(NA,6,12)
pos = which.max(CM.equil[,3])
msy.tab[1,] = c(Fv[pos],CM.boo.qRP[1,1],CM.boo.qRP[3,1],CM.boo.qRP[5,1],
                CM.equil[pos,2],CM.boo.qRP[1,2],CM.boo.qRP[3,2],CM.boo.qRP[5,2],
                CM.equil[pos,3],CM.boo.qRP[1,3],CM.boo.qRP[3,3],CM.boo.qRP[5,3]) 
pos = which.max(CM3.equil[,3])
msy.tab[2,] = c(Fv[pos],CM3.boo.qRP[1,1],CM3.boo.qRP[3,1],CM3.boo.qRP[5,1],
                CM3.equil[pos,2],CM3.boo.qRP[1,2],CM3.boo.qRP[3,2],CM3.boo.qRP[5,2],
                CM3.equil[pos,3],CM3.boo.qRP[1,3],CM3.boo.qRP[3,3],CM3.boo.qRP[5,3]) 
pos = which.max(CM4.equil[,3])
msy.tab[3,] = c(Fv[pos],CM4.boo.qRP[1,1],CM4.boo.qRP[3,1],CM4.boo.qRP[5,1],
                CM4.equil[pos,2],CM4.boo.qRP[1,2],CM4.boo.qRP[3,2],CM4.boo.qRP[5,2],
                CM4.equil[pos,3],CM4.boo.qRP[1,3],CM4.boo.qRP[3,3],CM4.boo.qRP[5,3])

msy.tab[4,] = c(hs.ref.points[3],HS.boo.qRP[1,3],HS.boo.qRP[3,3],HS.boo.qRP[5,3],
                hs.ref.points[4],HS.boo.qRP[1,4],HS.boo.qRP[3,4],HS.boo.qRP[5,4],
                hs.ref.points[5],HS.boo.qRP[1,5],HS.boo.qRP[3,5],HS.boo.qRP[5,5]) 
msy.tab[5,] = c(bh.ref.points[3],BH.boo.qRP[1,3],BH.boo.qRP[3,3],BH.boo.qRP[5,3],
                bh.ref.points[4],BH.boo.qRP[1,4],BH.boo.qRP[3,4],BH.boo.qRP[5,4],
                bh.ref.points[5],BH.boo.qRP[1,5],BH.boo.qRP[3,5],BH.boo.qRP[5,5]) 
msy.tab[6,] = c(rk.ref.points[3],RK.boo.qRP[1,3],RK.boo.qRP[3,3],RK.boo.qRP[5,3],
                rk.ref.points[4],RK.boo.qRP[1,4],RK.boo.qRP[3,4],RK.boo.qRP[5,4],
                rk.ref.points[5],RK.boo.qRP[1,5],RK.boo.qRP[3,5],RK.boo.qRP[5,5])
                

rownames(msy.tab) = c('CM GCV','CM df=3','CM df=4','HS','BH','RK')
colnames(msy.tab) = c('FMSY','LCI','Med','UCI','BMSY','LCI','Med','UCI','MSY','LCI','Med','UCI')

ctext = paste('MSY RP, Fmax = ',round(Fmax,digits=2),sep='')
print(xtable(msy.tab,digits=c(7,2,2,2,2,1,1,1,1,1,1,1,1),caption=ctext,
   align=c('c','r','r','r','r','r','r','r','r','r','r','r','r')),type='html',
  file='msy.html',caption.placement="top",hline.after=c(1,2))

##### Fit table #####

fit.tab = matrix(NA,6,3)                  
fit.tab[1,1] = sum((b.cm$weight*b.cm$residuals)**2)/n 
fit.tab[2,1] = sum((b.cm3$weight*b.cm3$residuals)**2)/n 
fit.tab[3,1] = sum((b.cm4$weight*b.cm4$residuals)**2)/n 
fit.tab[4:6,1] = apply((resid.matrix[,2:4])**2,2,mean) 

fit.tab[1,2] = sum(b.cm$edf)
fit.tab[2,2] = sum(b.cm3$edf)
fit.tab[3,2] = sum(b.cm4$edf)  
fit.tab[4,2] = 2        
fit.tab[5,2] = 2        
fit.tab[6,2] = 2

fit.tab[1,3] = mygcv(gcv.fit$par) 
fit.tab[2,3] = mygcv(log(sp3))        
fit.tab[3,3] = mygcv(log(sp4))  
fit.tab[4,3] = sum(HS.raw.res**2)/(n*(1 - 2/n)**2)        
fit.tab[5,3] = sum(BH.raw.res**2)/(n*(1 - 2/n)**2)        
fit.tab[6,3] = sum(RK.raw.res**2)/(n*(1 - 2/n)**2)

rownames(fit.tab) = c('CM GCV','CM df=3','CM df=4','HS','BH','RK')
colnames(fit.tab) = c('MSE','df','GCV')

ctext = paste('Fit statistics, n = ',n,sep='')
print(xtable(fit.tab,digits=c(7,3,1,3),caption=ctext,
  align=c('c','r','r','r')),type='html',
  file='fit_table.html',caption.placement="top")  
  
save.image('save.RData')

########################################################
