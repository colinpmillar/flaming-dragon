#
# this script simulates a stock recruitment relationship.
# Possibly coming from a rage of several models...
#
# Colin Millar January 10th 2013
#

# simulate from a gamma given the mean and cv
rgamma2 <- function(n, mean, cv) rgamma(n, 1/cv^2, 1/(mean*cv^2))


ssb <- seq(2, 10, length = 100)
par <- c(log(10), 0.2)
cv <- 0.3

X <- cbind(log(ssb), 1, -ssb)
invlink <- exp
ab <- c(1, par)

Erec <- invlink(X %*% ab)

rec <- rgamma2(length(ssb), Erec, cv)


plot(ssb, rec, xlim = c(0, max(ssb)), ylim = c(0, max(rec)))
lines(ssb, Erec)




