# Source models
source("models/model0.R")$value
source("models/model1.R")$value
source("models/model2.R")$value
source("models/model3.R")$value
source("models/model4.R")$value
source("models/model5.R")$value
source("models/model6.R")$value
source("models/model7.R")$value
source("models/model8.R")$value
source("models/model9.R")$value
# Source model parameters
source("parameters/param0.R")$value
source("parameters/param1.R")$value
source("parameters/param2.R")$value
source("parameters/param3.R")$value
source("parameters/param4.R")$value
source("parameters/param5.R")$value
source("parameters/param6.R")$value
source("parameters/param7.R")$value
source("parameters/param8.R")$value
source("parameters/param9.R")$value
# Source simulations
data0 <- read.csv("data/simulation0.csv", header = TRUE)
data1 <- read.csv("data/simulation1.csv", header = TRUE)
data2 <- read.csv("data/simulation2.csv", header = TRUE)
data3 <- read.csv("data/simulation3.csv", header = TRUE)
data4 <- read.csv("data/simulation4.csv", header = TRUE)
data5 <- read.csv("data/simulation5.csv", header = TRUE)
data6 <- read.csv("data/simulation6.csv", header = TRUE)
data7 <- read.csv("data/simulation7.csv", header = TRUE)
data8 <- read.csv("data/simulation8.csv", header = TRUE)
data9 <- read.csv("data/simulation9.csv", header = TRUE)
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


# Optimize model0
## Run optimizer
fit <- optim(par = param0, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data0[samps, ])),
             method = "Brent",
             lower = -10,
             upper = 10,
             disp.mod = model0,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimate values
t(param0) # true values


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


# Optimize model5
## Run optimizer
fit <- optim(par = param5, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data5[samps, ])),
             disp.mod = model5,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param5) # true values


# Optimize model6
## Run optimizer
fit <- optim(par = param6, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data6[samps, ])),
             disp.mod = model6,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param6) # true values


# Optimize model7
## Run optimizer
fit <- optim(par = param7, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data7[samps, ])),
             disp.mod = model7,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param7) # true values


# Optimize model8
## Run optimizer
fit <- optim(par = param8, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data8[samps, ])),
             disp.mod = model8,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param8) # true values


# Optimize model9
## Run optimizer
fit <- optim(par = param9, # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data9[samps, ])),
             disp.mod = model9,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param9) # true values

