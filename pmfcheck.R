# Source models
source("models/model1.R")$value
source("models/model2.R")$value
source("models/model3.R")$value
source("models/model4.R")$value
# Source model parameters
source("parameters/param1.R")$value
source("parameters/param2.R")$value
source("parameters/param3.R")$value
source("parameters/param4.R")$value
# Source simulations
data1 <- read.csv("data/simulation1.csv", header = TRUE)
data2 <- read.csv("data/simulation2.csv", header = TRUE)
data3 <- read.csv("data/simulation3.csv", header = TRUE)
data4 <- read.csv("data/simulation4.csv", header = TRUE)
# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently
# compute the exact pmf of a Poisson binomial distribution with function 'dpoibin()'
library(poibin)
# Source negative log-likelihood function (general)
source("nll.R")$value

# Randomly sample from simulations
## Set seed
set.seed(20190315)
## Set number of simulations
sims <- 100000
## Number of samples
num <- 500
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
fit$par # estimate values
t(param1) # true values

# Optimize model2
## Run optimizer
fit <- optim(par = param2, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data2[samps, ])),
             disp.mod = model2,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param2) # true values


# Optimize model3
## Run optimizer
fit <- optim(par = param3, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data3[samps, ])),
             disp.mod = model3,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param3) # true values


# Optimize model4
## Run optimizer
fit <- optim(par = param4, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data4[samps, ])),
             disp.mod = model4,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param4) # true values

