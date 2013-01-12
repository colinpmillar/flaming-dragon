######
# Main code for running MCMC algorithm for capture-recapture data
# - dippers when only considering model C/C. Created August 07 (RK)
# Uses single-update Uniform random walk step
######


######
# This function calculates the log likelihood of the stock recruit relationship
######

calcllikhood <- function(ni, data, param, model) {

# First we set the survival and recapture probs:

# Set up the size of the array containing the survival and recapture probs:
# Set up the set of cell probabilities (q), and initially set them to be all equal to zero:

    a <- param[1] 
    b <- param[2]
    cv <- param[3]

# Calculate the probability the log(likelihood):

    if (model == "bevholt") {
        prediction <- a * data $ ssb / (1 + b * data $ ssb)
    } else
    if (model == "ricker") {
        prediction <- a * data $ ssb * exp(-b * data $ ssb)
    } else
    if (model == "segreg") {
        prediction <- ifelse(data $ ssb < 1/(a*b), b * data $ ssb, 1/a)
    }
    llikhood <- sum( dnorm(log(data $ rec), log(prediction), cv, log = TRUE) )

# Output the log(likelihood) value:

    llikhood
}




######
# Propose a multiplicatively uniform variable ... (perhaps useful for MSE...)
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

updateparam <- function(nparam, param, ni, data, llikhood, alpha, beta, delta, model){

# Cycle through each parameter in turn and propose to update using
# random walk MH with Uniform proposal density:

    for (i in 1:nparam) {

# Update the parameter

# Keep a record of the current parameter value being updated

        oldparam <- param[i]

# Propose a new value for the parameter using a random walk with
# Uniform proposal density

        #param[i] <- param[1] + runif(1, -1, 1) * delta[i] # additive ... 
        param[i] <- param[1] * scaleProposal(delta[i])  # multiplicative

# Automatically reject any moves that are propose values less than zero

        if (param[i] > 0) {

# Calculate the log(acceptance probability):

# Calculate the new likelihood value for the proposed move:
# Calculate the numerator (num) and denominator (den) in turn:

            newllikhood <- calcllikhood(ni, data, param, model)

# Include the likelihood and prior (Beta) terms in the acceptance probability

            lnum <- newllikhood #+ dgamma(param[i],alpha[i],beta[i], log = TRUE)
            lden <- llikhood #+ dgamma(oldparam,alpha[i],beta[i], log = TRUE)

# All other prior terms (for other parameters) cancel in the acceptance probability.

# Acceptance probability of MH step:

            A <- min(1, exp(lnum - lden))

        } else {

            A <- 0

        }

# To do the accept/reject step of the algorithm:
# Simulate a random number in [0,1]:

        u <- runif(1)

# Accept the move with probability A:

        if (u <= A) {

# Accept the proposed move:
# Update the log(likelihood) value:

            llikhood <- newllikhood
        } else {

# Reject proposed move so parameter stays at current value:

            param[i] <- oldparam

        }
    }

# Set the values to be outputted from the function to be the
# set of parameter values and log(likelihood) value:

    output <- c(param, llikhood)

# Output the parameter values:

    output
}




updateparam <- function(nparam, param, ni, data, llikhood, alpha, beta, delta, model){

# Cycle through each parameter in turn and propose to update using
# random walk MH with Uniform proposal density:


# Update the parameter

# Keep a record of the current parameter value being updated

      oldparam <- param

# Propose a new value for the parameter using a random walk with
# Uniform proposal density

    #param[i] <- param[i] + runif(1, -1, 1) * delta[i] # additive ... 
    for (i in 1:nparam) {
        param[i] <- param[i] * scaleProposal(delta[1])  # multiplicative
    }

# Automatically reject any moves that are propose values less than zero

    if (all(param > 0)) {

# Calculate the log(acceptance probability):

# Calculate the new likelihood value for the proposed move:
# Calculate the numerator (num) and denominator (den) in turn:

        newllikhood <- calcllikhood(ni, data, param, model)

# Include the likelihood and prior (Beta) terms in the acceptance probability

        lnum <- newllikhood #+ dgamma(param[i],alpha[i],beta[i], log = TRUE)
        lden <- llikhood #+ dgamma(oldparam,alpha[i],beta[i], log = TRUE)

# All other prior terms (for other parameters) cancel in the acceptance probability.

# Acceptance probability of MH step:

        A <- min(1, exp(lnum - lden))
    } else {
        A <- 0
    }

# To do the accept/reject step of the algorithm:
# Simulate a random number in [0,1]:

    u <- runif(1)

# Accept the move with probability A:

    if (u <= A) {

# Accept the proposed move:
# Update the log(likelihood) value:

        llikhood <- newllikhood
    } else {

# Reject proposed move so parameter stays at current value:

        param <- oldparam

    }
  

# Set the values to be outputted from the function to be the
# set of parameter values and log(likelihood) value:

  output <- c(param, llikhood)

# Output the parameter values:

  output
}





# Define the function: with input parameters:
# nt = number of iterations
# nburn  = burn-in

BHMH <- function(nt, nburn, data, delta, model = c("bevholt","ricker","segreg")) {

# Define model

    model <- match.arg(model)

# Define the parameter values:
# ni = number of years
# nj = number of recapture years
# nparam = maximum number of parameters

    ni <- nrow(data)

# Set initial parameter values:

    if (model == "bevholt") {

      init <- glm(rec ~ I(1/ssb), family = Gamma(inverse), data = data)
      param. <- coef(init)
      param <- rep(NA,3)
      param[1] <- 1/param.[2]
      param[2] <- param.[1] * param[1]
      param[3] <- sqrt(gamma.dispersion(init))

    } else
    if (model == "ricker") {

      init <- glm(rec ~ ssb + offset(log(ssb)), family = Gamma(log), data = data)
      param. <- coef(init)
      param <- rep(NA,3)
      param[1] <- exp(param.[1])
      param[2] <- -1 * param.[2]
      param[3] <- sqrt(gamma.dispersion(init))

    } else
    if (model == "segreg") {

      func <- function(b, data) {
          init <- glm(rec ~ 1 + offset(ifelse(ssb < b, log(ssb / b), 0)), family = Gamma(log), data = data)
          -1 * unclass(logLik(init))
      }
      opt <- optimise(func, range(data $ ssb), data = data)
      b <- opt $ minimum
      init <- glm(rec ~ 1 + offset(ifelse(ssb < b, log(ssb / b), 0)), family = Gamma(log), data = data)

      param <- rep(NA,3)
      param[1] <- 1/exp(coef(init))
      param[2] <- exp(coef(init)) / b
      param[3] <- sqrt(gamma.dispersion( init ))
    }

    nparam <- length(param)

# param[1] = a
# param[2] = b
# param[3] = cv (sd on log scale...)

# sample is an array in which we put the sample from the posterior distribution.

    sample <- array(0, dim=c(nt, nparam + 1))
    sample <- data.frame(sample)
    names(sample) <- c("a", "b", "cv", "llik")
    accept <- rep(0, nparam)

# Calculate log(likelihood) for initial state using a separate function "calcllikhood":

    llikhood <- calcllikhood(ni, data, param, model)

# MCMC updates - MH algorithm:

# Cycle through each iteration:

    for (t in 1:nt) {

# Update the parameters in the model using function "updateparam":

        output <- updateparam(nparam, param, ni, data, llikhood, 
                              alpha, beta, delta, model)

# Set parameter values and log(likelihood) value of current state to be the output from
# the MH step:

        param <- output[1:nparam]
        llikhood <- output[nparam+1]

# Record the set of parameter values:

        sample[t,] <- output

    }

# Calculate the mean and standard deviation of the parameters
# following burn-in:

    subsample <- sample[(nburn+1):nt,]

    mn <- apply(subsample, 2, mean)
    std <- apply(subsample, 2, var)

# Output the posterior mean and standard deviation of the parameters
# following burn-in to the screen:

    # a crudish approximation
    acceptance.rate <- ( colSums((apply(subsample[1:3], 2, diff) != 0)) / nrow(subsample) )

    cat("Posterior summary estimates for each model:  ", "\n")
    cat("\n")
    cat("mean  (SD)", "\n")
    for (i in names(sample)[1:nparam]) {
        cat(i, ": ", "\n", sep = "")
        cat(mn[i], "   (", std[i], ") \tacc rate (try for 0.4) :", round(acceptance.rate[i], 3), "\n\n")
    }
# Output the sample from the posterior distribution:

    subsample
}


