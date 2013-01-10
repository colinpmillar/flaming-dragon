#*************************************************************************

#*************************************************************************
rm(list=ls())
# working directory

Anch <- c("C:/Documents and Settings/WinXP/Desktop/SGMED_2012_Ancona/Anchovy/ICA/Run5/")
setwd(Anch)

# species identifier

id <- c("A7511")

A11fl.dat <- readLines( paste(Anch,id,"fl.dat",sep="",collapse=""))
A11la.dat <- readLines( paste(Anch,id,"la.dat",sep="",collapse="") )

aux <- strsplit(A11la.dat[3], " ")[[1]]
aux <- aux[aux!=""]
aux <- as.numeric(aux)
minyear <- aux[1]
maxyear <- aux[2]
years <- minyear:maxyear
nyears <- length(years)
landings <- NULL
for (i in 1:nyears){
  aux <- strsplit(A11la.dat[5+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  landings <- c(landings, years[i], aux)
}
landings <- matrix(landings, ncol=2, byrow=T)

idx <- 1 # main title
idx <- 2 # number of surveys+100
aux <- strsplit(A11fl.dat[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- as.numeric(aux)
nageindices <- aux -100
ageindices <- NULL
for (k in 1:nageindices){
  idx <- idx + 2
  aux <- strsplit(A11fl.dat[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  minyear <- aux[1]
  maxyear <- aux[2]
  years <- minyear:maxyear
  nyears <- length(years)
  idx <- idx + 1
  aux <- strsplit(A11fl.dat[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  # aux[1] is sex code
  # aux[2] is effort code
  # aux[3] and aux[4] are start and end of survey as fractions of year
  idx <- idx + 1
  aux <- strsplit(A11fl.dat[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  minage <- aux[1]
  maxage <- aux[2]
  ages <- minage:maxage
  nages <- length(ages)
  nam <- paste(rep("Age", nages),c(minage:maxage),sep="")
  aux0 <- NULL
  for (i in 1:nyears){
    aux <- strsplit(A11fl.dat[idx+i], " ")[[1]]
    aux <- aux[aux!=""]
    aux <- as.numeric(aux)
    aux0 <- c(aux0, years[i], aux)
  }
  aux0[aux0==-1] <- NA
  aux0 <- matrix(aux0, nrow=nyears, byrow=T)
  aux0 <- aux0[,-2]
  aux0 <- as.data.frame(aux0)
  names(aux0) <- c("Year", nam)
  ageindices[[k]] <- aux0
  idx <- idx + nyears
}

nindices <- nageindices
#*************************************************************************
#*************************************************************************

# READ ICA FILES FROM PATTERSON

# read files for the output plots from ICA

# read icadiag.out
#
# contains estimates of the sum of squared residulas of the surveys compared with
# conventional VPA populations initiated with separable VPA, for 20 values of reference
# F over the range specified by the user 

icadiag.out <- readLines( paste(Anch,"icadiag.out",sep="",collapse="") )
nline <- length(icadiag.out)

ssq <- NULL
for (i in 4:nline){
  aux <- strsplit(icadiag.out[i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  ssq <- rbind(ssq, aux)
}
ssq <- as.data.frame(ssq, row.names=F)

# read ica.vie contains:
#
# contains:                                                  
# LANDINGS BY YEAR                                                                
# FISHING MORTALITY +/- STANDARD DEVIATION                                        
# RECRUITMENT                                                                     
# TOTAL AND SPAWNING BIOMASS                                                      
# SELECTION PATTERN +/- STANDARD DEVIATION                                        
# FISHING MORTALITY AT AGE AND YEAR                                               
# NATURAL MORTALITY BY AGE AND YEAR                                               
# BIOMASS INDEX CATCHABILITIES                                                    
# AGE-STRUCTURED INDEX CATCHABILITIES                                             
#

ica.vie <- readLines( paste(Anch,"/ica.vie",sep="",collapse="") )
nline <- length(ica.vie)

idx <- 2
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- as.numeric(aux)
minyear <- aux[1] # initial year in the analysis
maxyear <- aux[2] # final year in the analysis
nyears <- maxyear-minyear+1
minage <- aux[3] # initial age in the analysis
maxage <- aux[4] # last age in the analysis
nages <- maxage-minage+1
sepyear <- aux[5]
 # number of years in the separable analysis
# nssbindices <- aux[6] # number of ssb indices
# nageindices <- aux[7] # number of age-structured indices
# aux[8] no sé qué es.   

# landings are already read, so skip this part 

cont <- T
while(cont){
  idx <- idx +1
  aux <- strsplit(ica.vie[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
  cont <- (aux != "FISHING MORTALITY +/- STANDARD DEVIATION") & (idx < nline)
}

f.year <- NULL
for (i in 1:nyears){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  aux[aux<0] <- NA
  f.year <- rbind(f.year, aux)
}
row.names(f.year) <- NULL
f.year <- as.data.frame(f.year)
names(f.year) <- c("year","f","flow","fup")
idx <- idx + nyears + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="RECRUITMENT")
recruitment <- NULL
for (i in 1:(nyears+1)){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  recruitment <- rbind(recruitment, aux)
}
row.names(recruitment) <- NULL
recruitment <- as.data.frame(recruitment)
names(recruitment) <- c("year","r","rlow","rup")
idx <- idx + nyears + 1 + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="TOTAL AND SPAWNING BIOMASS")
biomass <- NULL
for (i in 1:(nyears+1)){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  biomass <- rbind(biomass, aux)
}
row.names(biomass) <- NULL
biomass <- as.data.frame(biomass)
names(biomass) <- c("year","tot","ssb")
idx <- idx + nyears + 1 + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="SELECTION PATTERN +/- STANDARD DEVIATION")
s.age <- NULL
for (i in 1:nages){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  aux[aux<0] <- NA
  s.age <- rbind(s.age, aux)
}

#s.age1 <- NULL
#for (i in 1:nages){
#  aux <- strsplit(ica.vie[idx+i+nages+1], " ")[[1]]
#  aux <- aux[aux!=""]
#  aux <- as.numeric(aux)
#  aux[aux<0] <- NA
#  s.age1 <- rbind(s.age1, aux)
#}
#
row.names(s.age) <- NULL
s.age <- as.data.frame(s.age)
names(s.age) <- c("age","s","slow","sup")
idx <- idx + nages + 1 #+ nages + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="NUMBERS AT AGE AND YEAR")
numbatage <- NULL
for (i in 1:(nyears+1)){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  numbatage <- rbind(numbatage, aux)
}
row.names(numbatage) <- NULL
idx <- idx + nyears + 1 + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="FISHING MORTALITY AT AGE AND YEAR")
f.ageyear <- NULL
for (i in 1:(nyears+1)){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  f.ageyear <- rbind(f.ageyear, aux)
}
row.names(f.ageyear) <- NULL

idx <- idx + nyears + 1 + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="NATURAL MORTALITY BY AGE AND YEAR")
m.ageyear <- NULL
for (i in 1:(nyears+1)){
  aux <- strsplit(ica.vie[idx+i], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  m.ageyear <- rbind(m.ageyear, aux)
}
row.names(m.ageyear) <- NULL

idx <- idx + nyears + 1 + 1
aux <- strsplit(ica.vie[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- paste(aux[aux!=""],sep=" ",collapse=" ")
print(aux =="BIOMASS INDEX CATCHABILITIES")

# FALTA LEER LAS CATCHABILITIES de ICAVIEW

icacn.res <- readLines( paste(Anch,"/icacn.res",sep="",collapse="") )
nline <- length(icacn.res)

idx <- 3
aux <- strsplit(icacn.res[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- as.numeric(aux)
minyear <- aux[1]
maxyear <- aux[2]

idx <- 4
aux <- strsplit(icacn.res[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- as.numeric(aux)
minage <- aux[1]
maxage <- aux[2]

catchatage.res <- NULL
for (idx in 6:nline){
  aux <- strsplit(icacn.res[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  catchatage.res <- rbind(catchatage.res, aux)
}
catchatage.res <- as.data.frame(catchatage.res, row.names=F)
catchatage.res <- as.matrix(catchatage.res)


# read file icaasurv.res

icaasurv.res <- readLines( paste(Anch,"/icaasurv.res",sep="",collapse="") )
nline <- length(icaasurv.res)

idx <- 1 # main title
idx <- 2 # number of surveys+100
aux <- strsplit(A11fl.dat[idx], " ")[[1]]
aux <- aux[aux!=""]
aux <- as.numeric(aux)

ageindices.res <- NULL
for (k in 1:nageindices){
  idx <- idx + 2
  aux <- strsplit(icaasurv.res[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  minyear <- aux[1]
  maxyear <- aux[2]
  years <- minyear:maxyear
  nyears <- length(years)
  idx <- idx + 1
  aux <- strsplit(icaasurv.res[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  # aux[1] is sex code
  # aux[2] is effort code
  # aux[3] and aux[4] are start and end of survey as fractions of year
  idx <- idx + 1
  aux <- strsplit(icaasurv.res[idx], " ")[[1]]
  aux <- aux[aux!=""]
  aux <- as.numeric(aux)
  minage <- aux[1]
  maxage <- aux[2]
  ages <- minage:maxage
  nages <- length(ages)
  aux0 <- NULL
  for (i in 1:nyears){
    aux <- strsplit(icaasurv.res[idx+i], " ")[[1]]
    aux <- aux[aux!=""]
    aux <- as.numeric(aux)
    aux[aux==99.99] <- NA
    aux0 <- c(aux0, aux)
  }
  aux0[aux0==-1] <- NA
  aux0 <- matrix(aux0, nrow=nyears, byrow=T)
  aux0 <- as.data.frame(aux0)
  ageindices.res[[k]] <- aux0
  idx <- idx + nyears
}

rm(icadiag.out, ica.vie, icacn.res, icabsurv.res, icaasurv.res)

#*************************************************************************

# open pdf to write out plots

pdf(paste(Anch, id, "_run5.pdf", sep="",collapse=""))

# fig1

matplot(ssq[,1], ssq[, -1], type="l", xlab="Reference F", ylab="Index SSQ", main="SSQ Surface", col=nindices, lty=nindices)
legend(0, max(ssq[, -1]), paste("Age"), bty="n", col=nindices, lty=nindices, cex=1, xjust=-0.5)

# fig2

par(mfrow=c(2,2))

plot(landings, type="l", xlab="Year", ylab="Yield", main="Landings")

plot(range(f.year[1:dim(f.year)[1],1]), range(f.year[1:dim(f.year)[1], -1], na.rm=T), type="n", xlab="Year", ylab="Reference F (1-3)", yaxp=c(0,2,10), main="Fishing mortality")
lines(f.year[1:dim(f.year)[1],1:2])
for (i in 1:dim(f.year)[1]){
  segments(f.year[i,1], f.year[i,3], f.year[i,1], f.year[i,4],col="darkblue", lty=2)
  segments(f.year[i,1]-0.2, f.year[i,3], f.year[i,1]+0.2, f.year[i,3], col="darkblue", lty=2)
  segments(f.year[i,1]-0.2, f.year[i,4], f.year[i,1]+0.2, f.year[i,4],col="darkblue", lty=2)
}

plot(recruitment[1:dim(recruitment)[1], 1:2], type="l", xlim=c(1975,2011), xlab="Year", ylab="Recruits", main="Recruitment")

plot(range(biomass[1:dim(biomass)[1],1]), range(biomass[1:dim(biomass)[1], 3], na.rm=T), type="n", xlab="Year", ylab="Biomass", main="Stock size (Mid Year)")
#lines(biomass[1:dim(biomass)[1],1], biomass[1:dim(biomass)[1],2])
lines(biomass[1:dim(biomass)[1],1], biomass[1:dim(biomass)[1],3], lty=1)
#legend(min(biomass[1:dim(biomass)[1],1]), max(biomass[1:dim(biomass)[1],-1], na.rm=T), c("Total","SSB"), lty=1:2, cex=0.6)

# fig3

par(mfrow=c(2,2))

contour(2000:(2000+dim(catchatage.res)[1]-1), 0:(dim(catchatage.res)[2]-1), catchatage.res, xlab="Year", ylab="Age", main="Log Residual")
#contour(minyear:maxyear, minage:maxage, catchatage.res, xlab="Year", ylab="Age", main="Log Residual")

plot(range(s.age[,1]), range(s.age[, -1], na.rm=T), type="n", xlab="Age", ylab="Selection", main="Selection pattern")
lines(s.age[,1:2])
lines(s.age1[,1:2])
for (i in 1:dim(s.age)[1]){
  segments(s.age[i,1], s.age[i,3], s.age[i,1], s.age[i,4])
  segments(s.age[i,1]-0.2, s.age[i,3], s.age[i,1]+0.2, s.age[i,3])
  segments(s.age[i,1]-0.2, s.age[i,4], s.age[i,1]+0.2, s.age[i,4])
}


barplot(apply(catchatage.res, 1, sum), names.arg=2002:2011, space=0, cex.names=0.7, ylim=c(-1,+1), xlab="Year", ylab="Marginal total", main="Year Residuals")

barplot(apply(catchatage.res, 2, sum), names.arg=0:5, space=0, cex.names=0.7, ylim= c(-1, +1), xlab="Age", ylab="Marginal total", main="Age Residuals")


#fig5

#par(mfrow=c(2,4))

#for (j in 1:nageindices){
#  for (k in 1: (dim(ageindices[[j]])[2]-1) ){
#    plot(ageindices[[j]][,1], ageindices.res[[j]][,k], xlab="Year", ylab="Residual")
#    abline(h=0, lty=3)
#  }
#}

dev.off()

#*************************************************************************
#*************************************************************************


