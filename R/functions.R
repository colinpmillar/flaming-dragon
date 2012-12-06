

BH <- function(parm,s){
  alpha <- parm[1]
  beta <- exp(parm[2]);
  ksi = alpha + log(s) - log(beta + s);
  return(exp(ksi))
}

BH.ll <- function(parm){
  s<-x
  ly <- log(y)
  res=sum(wt.ext*((ly-log(BH(parm,s)))^2))
  return(res)
}


BH.gam <- function(parm){
  mu = BH(parm,x);
  temp = y/mu
  res=sum(2*(temp - log(temp) -1))
  return(res)
}

RK <- function(parm,s){
  alpha <- parm[1]
  beta <- exp(parm[2]);
  ksi = alpha + log(s) - beta*s;
  return(exp(ksi))
}

RK.ll <- function(parm) {
  s<-x
  ly <- log(y)
  res=sum(wt.ext*((ly-log(RK(parm,s)))^2))
  return(res)
}


RK.gam <- function(parm){
  mu = RK(parm,x);
  temp = y/mu
  res=sum(2*(temp - log(temp) -1))
  return(res)
}


fS50 <- function(sp){
  return((hmax.rec - RK(tparm,sp))^2)
}


HS <- function(parm,s){
  ls = log(s)
  alpha <- parm[1]
  delta <- parm[2];
  log.rec = alpha + log(s)
  ind = ls>delta
  log.rec[ind] = alpha + delta
  return(exp(log.rec))
}


HS.ll <- function(parm){
  s<-x
  ly <- log(y)
  res=sum(wt.ext*((ly-log(HS(parm,s)))^2))
  return(res)
}


HS.gam <- function(parm){
  mu = HS(parm,x);
  temp = y/mu
  res=sum(2*(temp - log(temp) -1))
  return(res)
}

se <- function(H){
  det.h <- 1/det(H)
  if(is.na(det.h)){se <- NA}
  if(!is.na(det.h)){
  if(det.h>1e-3){
    cma <- chol(H)
    se <- sqrt(diag(chol2inv(cma)))
  }
  if(det.h<1e-3){
    se <- c(sqrt(1/H[1,1]),NA)
  }}
return(se)
}



mixBH <- function(parm,s){
  alpha <- parm[1]
  beta <- exp(parm[2]);
  mu1 = exp(alpha + log(s) - log(beta + s));
  alpha <- parm[3]
  beta <- exp(parm[4]);
  mu2 = exp(alpha + log(s) - log(beta + s));
  ## An mix range (75%-25%) of 200 ssb unit;
  logit = log(9)*(s - exp(parm[5]))/200
  mix = exp(logit)/(1+exp(logit))
  mu = mix*mu1 + (1-mix)*mu2;
  ret = cbind(mu1,mu2,mu,mix)
  return(ret)
}


mixBH.gam <- function(parm){
  mu = mixBH(parm,x)[,3];
  temp = y/mu
  res=sum(2*(temp - log(temp) -1))
  return(res)
}


## MSY Stuff;

equil.res <- function(f, alpha.p, beta.p, Saop, S50p, Rmaxp, sr_type) {
    fp = f * pr
    zp = mp + fp
    czp = cumsum( c(0, zp[1:(n.age-1)]))

    YPR = sum(exp(-czp) * catch.weight * fp * (1 - exp(-zp)) / zp)
#    if(!mid.year) {
      BPR = sum(exp(-czp)*stock.weight**matp)
#    } else {
#      BPR = sum(exp(-czp - zp/2)*stock.weight*matp)} ## midyear SSB per recruit
#    }
    BPRi = 1/BPR
    if(sr_type == 'BH') {B.equil = S50p * (Saop * BPR - 1)}
    if(sr_type == 'RK') {B.equil = (log(alpha.p) + log(BPR)) / beta.p}
    if(sr_type == 'HS') {
      B.equil = 0
      if (BPRi < Saop) {B.equil = Rmaxp/BPRi}
    }
    B.equil = max(0, B.equil)
    R.equil = B.equil / BPR
    Y.equil = R.equil * YPR

  c(R.equil, B.equil, Y.equil)
}

equil.yield <- function(f, alpha.p, beta.p, Saop, S50p, Rmaxp, sr_type) {
  -equil.res(f, alpha.p, beta.p, Saop, S50p, Rmaxp, sr_type)[3]
}

fS50 <- function(sp,hmax.rec,rparm){
  (hmax.rec - RK(rparm,sp))^2
}

make.RP <- function(parm, Fstart, iter.max = 200, sr_type){
  alpha.p = exp(parm[1])
  beta.p = exp(parm[2])
  upper = 4
  if(sr_type == 'BH'){
    S50p = exp(parm[2])
    Saop = exp(parm[1]-parm[2])
    Rmaxp = exp(parm[1])
  }
  if(sr_type == 'RK'){
   Saop = exp(parm[1])
   Rmaxp = alpha.p/(beta.p*exp(1))
   hmax.rec = Rmaxp/2
   ts50 = 1/(4.44*beta.p)
   rparm = parm;
   temp1 = nlminb(ts50,fS50,hmax.rec=hmax.rec,rparm=rparm)
   S50p = temp1$par
  }
  if(sr_type == 'HS'){
    S50p = beta.p/2
    Saop = 2*alpha.p
    Rmaxp = alpha.p*(beta.p + sqrt(beta.p**2 + gam2by4))
    fv = seq(0.2,2,by=0.001)
    mfv = matrix(fv,ncol=1)
    yc = apply(mfv,1,equil.yield,
         alpha.p=alpha.p,beta.p=beta.p,Saop=Saop,S50p=S50p,Rmaxp=Rmaxp, sr_type=sr_type);
    pos = which.min(yc)
    if(iter.max==0){upper = 4}
    if(iter.max>0){upper = mfv[pos,1]}
  }
  #temp = nlminb(Fstart,equil.yield, control = list(trace=1))
  temp = nlminb(Fstart, equil.yield, control = list(iter.max = iter.max), upper = upper,
         alpha.p = alpha.p, beta.p = beta.p, Saop = Saop, S50p = S50p, Rmaxp = Rmaxp, sr_type = sr_type)
  temp1 = equil.res(temp$par, alpha.p, beta.p, Saop, S50p, Rmaxp, sr_type)
  c(S50p,Saop,temp$par,temp1[2:3])
}


