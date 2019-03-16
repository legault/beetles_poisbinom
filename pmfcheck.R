# Source models
source("models/model1.R")
source("models/model2.R")
# Source model parameters
source("parameters/param1.R")
source("parameters/param2.R")
# Source simulations
data1 <- read.csv("data/simulation1.csv", header = TRUE)
data2 <- read.csv("data/simulation2.csv", header = TRUE)
# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently
# compute the exact pmf of a Poisson binomial distribution with function 'dpoibin()'
library(poibin)
# Source negative log-likelihood function (general)
source("nll.R")

# Randomly sample from simulations
## Set seed
set.seed(20190315)
## Set number of simulations
sims <- 100000
## Number of samples
num <- 200
## Generate samples
samps <- sample(1:sims, size = num, replace = TRUE)

# Optimize model1
## Run optimizer
fit <- optim(par = param1, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data1[samps, ])),
             disp.mod = model1,
             control = list(trace = 1))
## Summarize fit
fit
fit$par
## True values
t(param1)

# Optimize model2
## Run optimizer
fit <- optim(par = param2, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data2[samps, ])),
             disp.mod = model2,
             control = list(trace = 1))
## Summarize fit
fit
fit$par
## True values
t(param2)
