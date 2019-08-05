# Source models
sapply(list.files(path = "models", pattern = "*.R"), FUN = function(X) source(paste("models/", X, sep = "")))
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
data10 <- read.csv("data/simulation10.csv", header = TRUE)
data11 <- read.csv("data/simulation11.csv", header = TRUE)
data12 <- read.csv("data/simulation12.csv", header = TRUE)
data13 <- read.csv("data/simulation13.csv", header = TRUE)
data14 <- read.csv("data/simulation14.csv", header = TRUE)
data15 <- read.csv("data/simulation15.csv", header = TRUE)
# Source parameters
source("simulations/param.R")

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
num <- 500
## Generate samples
samps <- sample(1:sims, size = num, replace = TRUE)

# Optimize model0
## Run optimizer
fit <- optim(par = param[1], # use true values as initial values
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
t(param[1]) # true values

# Optimize model1
## Run optimizer
fit <- optim(par = param[1:2], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data1[samps, ])),
             disp.mod = model1,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimate values
t(param[1:2]) # true values

# Optimize model2
## Run optimizer
fit <- optim(par = param[c(1:3)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data2[samps, ])),
             disp.mod = model2,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1:3)]) # true values

# Optimize model3
## Run optimizer
fit <- optim(par = param[c(1, 2, 4)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data3[samps, ])),
             disp.mod = model3,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 2, 4)]) # true values

# Optimize model4
## Run optimizer
fit <- optim(par = param[c(1, 2, 4, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data4[samps, ])),
             disp.mod = model4,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 2, 4, 5)]) # true values

# Optimize model5
## Run optimizer
fit <- optim(par = param[1:5], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data5[samps, ])),
             disp.mod = model5,
             control = list(trace = 1, maxit = 1000)) # requires a few more iterations
## Summarize fit
fit
fit$par # estimated values
t(param[1:5]) # true values

# Optimize model6
## Run optimizer
fit <- optim(par = param[c(1, 4)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data6[samps, ])),
             disp.mod = model6,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 4)]) # true values

# Optimize model7
## Run optimizer
fit <- optim(par = param[c(1, 4, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data7[samps, ])),
             disp.mod = model7,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 4, 5)]) # true values

# Optimize model8
## Run optimizer
fit <- optim(par = param[c(1, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data8[samps, ])),
             disp.mod = model8,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 5)]) # true values

# Optimize model9
## Run optimizer
fit <- optim(par = param[c(1, 2, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data9[samps, ])),
             disp.mod = model9,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 2, 5)]) # true values

# Optimize model10
## Run optimizer
fit <- optim(par = param[c(1, 3)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data10[samps, ])),
             disp.mod = model10,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 3)]) # true values

# Optimize model11
## Run optimizer
fit <- optim(par = param[c(1, 3, 4)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data11[samps, ])),
             disp.mod = model11,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 3, 4)]) # true values

# Optimize model12
## Run optimizer
fit <- optim(par = param[c(1:4)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data12[samps, ])),
             disp.mod = model12,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1:4)]) # true values

# Optimize model13
## Run optimizer
fit <- optim(par = param[c(1, 3, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data13[samps, ])),
             disp.mod = model13,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 3, 5)]) # true values

# Optimize model14
## Run optimizer
fit <- optim(par = param[c(1, 2, 3, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data14[samps, ])),
             disp.mod = model14,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 2, 3, 5)]) # true values

# Optimize model15
## Run optimizer
fit <- optim(par = param[c(1, 3, 4, 5)], # use true values as initial values
             fn = nll,
             data = as.list(data.frame(data15[samps, ])),
             disp.mod = model15,
             control = list(trace = 1))
## Summarize fit
fit
fit$par # estimated values
t(param[c(1, 3, 4, 5)]) # true values

