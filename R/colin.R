
library(colorout)
library(FLa4av2)


stk <- readFLStock("S7511index.dat")
catch.n(stk) <- readVPAFile("S7511cn.dat")
catch.wt(stk) <- readVPAFile("S7511cw.dat")

#stk <- window(stk, start = 1980, end = 2011)

ind <- FLCore:::readIndicesVPA("S7511flno0.dat", quiet = FALSE, sep = "", na.strings = NA)
ind <- window(ind[[1]], start = 2004, end = 2011)

index(ind)[index(ind) == 0] <- 2000

#fmodel <- ~ s(year, k = 7, by = age)
fmodel <- ~ s(age, k = 3) + factor(year)

qmodel <- list(~ s(age, k=3))
#extra <- list(age1 = expression(as.numeric(age==0)),
#              age2 = expression(as.numeric(age==2)))

rmodel <- ~ factor(year)

fit <- a4aFit(fmodel, qmodel, rmodel, stk, list(ind))



