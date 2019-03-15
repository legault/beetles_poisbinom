# Parameters
sims <- 100000 # number of simulations
p0 <- log(.1 / (1 - .1)) # density-independent probability of dispersing (logit scale)
b1 <- .01 # density-dependent increase in probability of dispersing (logit scale)

# Source model1
source("models/model1.R")

# Simulate dispersal
## Set seed
set.seed(20190314)
## Create storage matrix
sim.stor <- matrix(NA, ncol = 4, nrow = sims)
## Fill matrix
sim.stor[, 1] <- rpois(n = sims, lambda = 80) # random abundance before dispersal in patch1
sim.stor[, 2] <- rpois(n = sims, lambda = 50) # random abundance before dispersal in patch2
d1 <- rbinom(n = sims, size = sim.stor[, 1], prob = model1(sim.stor[, 1], p0, b1)) # dispersal patch1
d2 <- rbinom(n = sims, size = sim.stor[, 2], prob = model1(sim.stor[, 2], p0, b1)) # dispersal patch2
sim.stor[, 3] <- sim.stor[, 1] - d1 + d2 # abundance after dispersal in patch1
sim.stor[, 4] <- sim.stor[, 2] - d2 + d1 # abundance after dispersal in patch2
## Rename columns
colnames(sim.stor) <- c("patch1.before", "patch2.before", "patch1.after", "patch2.after")
## Check sims
head(sim.stor)


# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently compute the exact pmf of a Poisson binomial distribution (function 'dpoibin()')
library(poibin)
    
nll <- function(param, data, disp.mod){
    d1 <- do.call(disp.mod, append(data, param)) # probability disperse from patch1
    d2 <- do.call(disp.mod, append(data, param)) # probability disperse from patch2
    loglike <- dim(length(data$patch1.before))
    for(i in 1:length(data$patch1.before)){
        loglike[i] <- log(dpoibin(kk = data$patch1.after[i],
                                  pp = c(rep(1 - d1[i], data$patch1.before[i]),
                                         rep(d2[i], data$patch2.before[i]))))
    }
    loglike[is.na(loglike)] <- -10000 # punish Inf
    return(-sum(loglike))
}

# Randomly sample from simulations
## Set seed
set.seed(20190315)
## Number of samples
num <- 200
## Generate samples
samps <- sample(1:sims, size = num, replace = TRUE)

# Optimize
## Set starting values for optimization
init <- list(p0 = p0,
             b1 = b1)
## Run optimizer
fit <- optim(par = init,
             fn = nll,
             data = as.list(data.frame(sim.stor[samps, ])),
             disp.mod = model1,
             control = list(trace = 1))
## Summarize fit
fit
p0
b1
fit$par
