


##############################################
# stock recruit formulations
##############################################


ricker <- function(ab, ssb) {
  log(ab$a) + log(ssb) - ab$b * ssb
}

segreg <- function(ab, ssb) {
  log(ifelse(ssb >= ab$b, ab$a * ab$b, ab$a * ssb))
}

bevholt <- function(ab, ssb) {
  log(ab$a * ssb / (1 + ab$b * ssb))
}


segreg2 <- function(ab, ssb) {
  log(ifelse(ssb >= ab$b, ab$a, ab$a / ab$b * ssb))
}

bevholt2 <- function(ab, ssb) {
  log(ab$a * ssb / (ab$b + ssb))
}



######
# This function calculates the log likelihood of the stock recruit relationship
######

llik <- function(param, data, model, logpar = FALSE) 
{
  FUN <- match.fun(model)
  if (logpar) {
    pred <- FUN(list(a = exp(param[1]), b = exp(param[2])), data $ ssb)
    sum( dnorm(log(data $ rec), pred, exp(param[3]), log = TRUE) )
  } else {
    pred <- FUN(list(a = param[1], b = param[2]), data $ ssb)
    sum( dnorm(log(data $ rec), pred, param[3], log = TRUE) )
  }
}


######
# Propose a multiplicatively uniform variable ...
######
scaleProposal <- function(A)
{
    if (A <= 1.0) {
        1.0
    } else {
        len = A - 1 / A
        if (runif(1) < len / (len + 2 * log(A))) {
            1 / A + len * runif(1)
        } else {
            A^(2.0 * runif(1) - 1)
        }
    }
}


######
# Function for updating the parameters values:
######

updateparam <- function(param, data, llikhood, delta, model)
{
# Keep a record of the current parameter value being updated

  oldparam <- param
  for (i in seq(param)) {
    param[i] <- param[i] * scaleProposal(delta)  # multiplicative
  }  
  newllikhood <- llik(param, data, model)

# MH step:

  if (runif(1) <= exp(newllikhood - llikhood)) { # no need for min(0, ...) 
# Accept the proposed move:
    llikhood <- newllikhood
  } else {
    param <- oldparam
  }

  c(llikhood, param)
}


MH <- function(nt, nburn, data, model, delta = 1.3) {

# scale data
  sdata <- data
  sdata $ ssb <- sdata $ ssb / exp(mean(log(sdata $ ssb)))
  sdata $ rec <- sdata $ rec / exp(mean(log(sdata $ rec)))

# Set initial parameter values - use maximum likelihood estimates:

  opt <- optim(rep(0, 3), llik, data = sdata, model = model, logpar = TRUE, control = list(fnscale = -1)) 
  param <- exp(opt $ par)

# param[1] = a
# param[2] = b
# param[3] = cv (sd on log scale...)

# sample is an array in which we put the sample from the posterior distribution.

  sample <- array(0, dim=c(nt, 4))
  sample <- data.frame(sample)
  names(sample) <- c("llik", "a", "b", "cv")

# Calculate log(likelihood) for initial state using a separate function "llik":

  llikhood <- llik(param, sdata, model)

# MCMC updates - MH algorithm:

  for (t in 1:nt) 
  {
    output <- updateparam(param, sdata, llikhood, delta, model)

# Set parameter values and log(likelihood) value of current state to be the input to
# the next MH step:
    param <- output[-1]
    llikhood <- output[1]

# Record the set of parameter values:
    sample[t,] <- output
  }

# Calculate the mean and standard deviation of the parameters
# following burn-in:
  subsample <- sample[(nburn+1):nt,]

# a crudish approximation to acceptance rate
  cat("acceptance rate:", round(mean(diff(subsample[,1]) != 0), 3), ", try for 0.40\n")

  subsample
}


