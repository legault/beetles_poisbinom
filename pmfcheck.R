# Parameters
sims <- 10000 # number of simulations
p0 <- log(.1 / (1 - .1)) # density-independent probability of dispersing (logit scale)
b1 <- .01 # density-dependent increase in probability of dispersing (logit scale)

pdisp <- function(p0, b1, x){
    lmod <- p0 + b1 * x
    return(exp(lmod) / (exp(lmod) + 1))
}

# Simulate dispersal
## Set seed
set.seed(20190314)
## Create storage matrix
sim.stor <- matrix(NA, ncol = 4, nrow = sims)
## Fill matrix
sim.stor[, 1] <- rpois(n = sims, lambda = 80) # random abundance before dispersal in patch1
sim.stor[, 2] <- rpois(n = sims, lambda = 50) # random abundance before dispersal in patch2
d1 <- sapply(1:sims, FUN = function(X){rbinom(n = 1, size = sim.stor[X, 1], prob = pdisp(p0, b1, sim.stor[X, 1]))}) # dispersal from patch1
d2 <- sapply(1:sims, FUN = function(X){rbinom(n = 1, size = sim.stor[X, 2], prob = pdisp(p0, b1, sim.stor[X, 2]))})
sim.stor[, 3] <- sim.stor[, 1] - d1 + d2 # abundance after dispersal in patch1
sim.stor[, 4] <- sim.stor[, 2] - d2 + d1 # abundance after dispersal in patch2
## Rename columns
colnames(sim.stor) <- c("patch1.before", "patch2.before", "patch1.after", "patch2.after")
## Check sims
head(sim.stor)


# Source 'poibin' library, which uses the discrete fourier tranform (DFT) to exactly compute the cdf and pmf of a Poisson binomial distribution
library(poibin)
    
nll <- function(param, patch1.before, patch2.before, patch1.after){
    p0 <- param[[1]]
    b1 <- param[[2]]
    d1 <- pdisp(p0, b1, patch1.before) # probability of dispersing from patch1
    d2 <- pdisp(p0, b1, patch2.before) # probability of dispersing from patch2
    loglike <- dim(length(patch1.before))
    for(i in 1:length(patch1.before)){
        loglike[i] <- log(dpoibin(kk = patch1.after[i],
                                  pp = c(rep(d2[i], patch2.before[i]),
                                         rep(1 - d1[i], patch1.before[i]))))
    }
    loglike[is.na(loglike)] <- -10000 # punish Inf
    return(-sum(loglike))
}

                        
init <- c(p0 = p0,
          b1 = b1)

# Randomly sample from simulations
num <- 800 # number of samples
samps <- sample(1:sims, size = num, replace = TRUE)
# Optimize
but <- optim(par = init,
             fn = nll,
             patch1.before = sim.stor[samps, 1],
             patch2.before = sim.stor[samps, 2],
             patch1.after = sim.stor[samps, 3],
             control = list(trace = 1))
but
p0
b1
but$par
    
## pmf.test <- function(k, n1, n2, p1, p2){
##     C <- exp((2 * 1i * pi) / (n1 + n2 + 1))
##     innerP <- prod(sapply(1:n1, FUN = function(m){1 + (C ^ m - 1) * p1}),
##                    sapply(1:n2, FUN = function(m){1 + (C ^ m - 1) * p2}))
##     outerS <- sum(sapply(0:(n1 + n2), FUN = function(l){C ^ (-l * k) * innerP}))
##     (1 / (n1 + n2 + 1)) * outerS
## }

## }


## pmf.test <- function(k, n1, n2, p1, p2){
##     omega <- (2 * pi) / (n1 + n2 + 1)
##     z1 <- function(l, p1){
##         1 - p1 + (p1 * cos(omega * l)) + (1i * p1 * sin(omega * l))
##     }
##     z2 <- function(l, p2){
##         1 - p2 + (p2 * cos(omega * l)) + (1i * p2 * sin(omega * l))
##     }
## }


## Arg <- function(l, p, omega){
##         atan2(y = p * sin(omega * l),
##               x = 1 - p + (p * cos(omega * l)))
## }

## modz <- function(l, p, omega){
##     sqrt((1 - p + (p * cos(omega * l))) ^ 2 + (p * sin(omega * l)) ^ 2)
## }


## pmf.test <- function(k, n1, n2, p1, p2){
##     omega <- (2 * pi) / (n1 + n2 + 1)
##     d.l <- sum(sapply(
    
        
