# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently
# compute the exact pmf of a Poisson binomial distribution with function 'dpoibin()'
library(poibin)
# Source plyr for rbind.full function
library(plyr)
# Source xtable for saving tables
library(xtable)

# Create list of models
m.list <- list.files(path = "models", pattern = "*.R")
# Source pmf functions
sapply(m.list, FUN = function(X) source(paste("models/", X, sep = "")))
# Trim ".R" extension from list
m.list <- gsub(pattern = "\\.R$", "", m.list)

# Set seed
set.seed(20190422)

# Source negative log-likelihood function (general)
source("nll.R")$value

# Source (rounded) data
csdata <- read.csv("data/csdata.csv", header = TRUE)
cfdata <- read.csv("data/cfdata.csv", header = TRUE)

# Create empty table of estimates
estimates <- data.frame(Species = c(),
                        Model = c(),
                        Parameters = c(),
                        negloglike = c(),
                        set = c(),
                        converge = c(),
                        p0 = c(),
                        b1 = c(),
                        b2 = c(),
                        b3 = c(),
                        b4 = c())

# CS
## Constant dispersal rate
startingvalues <- list(p0 = -1)
### Run optimizer
for(i in 1){
    arg.temp <- formalArgs(paste(m.list[i]))
    startingvalues.temp <- startingvalues[which(names(startingvalues) %in% arg.temp)]
    fit <- optim(par = startingvalues.temp,
                 fn = nll,
                 data = as.list(csdata),
                 disp.mod = paste(m.list[i]),
                 method = c("Brent"),
                 lower = -10,
                 upper = 10,
                 control = list(maxit = 1000, trace = 1))
    estimates.temp <- cbind(data.frame(Species = c("CS"),
                                       Model = m.list[i],
                                       Parameters = length(fit$par),
                                       negloglike = fit$value,
                                       set = 1,
                                       converge = fit$convergence,
                                       p0 = fit$par))
    estimates <- rbind.fill(estimates, estimates.temp)
    startingvalues.temp <- as.list(fit$par)
}

## All other dispersal models
### Number of random starting values
numr <- 200
### Use estimated p0 in subsequent optimization; random values for other parameters
startingvalues <- list(p0 = rep(fit$par, numr),
                       b1 = rnorm(n = numr, mean = 0, sd = .01),
                       b2 = rnorm(n = numr, mean = 0, sd = .01),
                       b3 = rnorm(n = numr, mean = 0, sd = .01),
                       b4 = rnorm(n = numr, mean = 0, sd = .01))
### Run optimizer
for(i in 2:length(m.list)){
    arg.temp <- formalArgs(paste(m.list[i]))
    startingvalues.i <- startingvalues[which(names(startingvalues) %in% arg.temp)]
    for(j in 1:numr){
        startingvalues.j <- lapply(startingvalues.i, "[[", j)
        fit <- optim(par = startingvalues.j,
                     fn = nll,
                     data = as.list(csdata),
                     disp.mod = paste(m.list[i]),
                     method = c("Nelder-Mead"),
                     control = list(maxit = 2000))
        if(fit$convergence == 10){
            t <- 1
            maxt <- 10
            while(t <= maxt){
                fit <- optim(par = lapply(fit$par, FUN = function(X) rnorm(n=1, mean=X, sd=0.01)),
                             fn = nll,
                             data = as.list(csdata),
                             disp.mod = paste(m.list[i]),
                             method = c("Nelder-Mead"),
                             control = list(maxit = 2000))
                if(fit$convergence == 0) break
                t <- t + 1
            }
        }
        estimates.temp <- cbind(data.frame(Species = c("CS"),
                                           Model = m.list[i],
                                           Parameters = length(fit$par),
                                           negloglike = fit$value,
                                           set = j,
                                           converge = fit$convergence),
                                t(as.data.frame(fit$par)))
        estimates <- rbind.fill(estimates, estimates.temp)
        startingvalues.temp <- as.list(fit$par)
    }
}

# CF
## Constant dispersal rate
startingvalues <- list(p0 = -1)
### Run optimizer
for(i in 1){
    arg.temp <- formalArgs(paste(m.list[i]))
    startingvalues.temp <- startingvalues[which(names(startingvalues) %in% arg.temp)]
    fit <- optim(par = startingvalues.temp,
                 fn = nll,
                 data = as.list(cfdata),
                 disp.mod = paste(m.list[i]),
                 method = c("Brent"),
                 lower = -10,
                 upper = 10,
                 control = list(maxit = 1000, trace = 1))
    estimates.temp <- cbind(data.frame(Species = c("CF"),
                                       Model = m.list[i],
                                       Parameters = length(fit$par),
                                       negloglike = fit$value,
                                       set = 1,
                                       converge = fit$convergence,
                                       p0 = fit$par))
    estimates <- rbind.fill(estimates, estimates.temp)
    startingvalues.temp <- as.list(fit$par)
}

## All other dispersal models
### Number of random starting values
numr <- 200
### Use estimated p0 in subsequent optimization; random values for other parameters
startingvalues <- list(p0 = rep(fit$par, numr),
                       b1 = rnorm(n = numr, mean = 0, sd = .01),
                       b2 = rnorm(n = numr, mean = 0, sd = .01),
                       b3 = rnorm(n = numr, mean = 0, sd = .01),
                       b4 = rnorm(n = numr, mean = 0, sd = .01))
### Run optimizer
for(i in 2:length(m.list)){
    arg.temp <- formalArgs(paste(m.list[i]))
    startingvalues.i <- startingvalues[which(names(startingvalues) %in% arg.temp)]
    for(j in 1:numr){
        startingvalues.j <- lapply(startingvalues.i, "[[", j)
        fit <- optim(par = startingvalues.j,
                     fn = nll,
                     data = as.list(cfdata),
                     disp.mod = paste(m.list[i]),
                     method = c("Nelder-Mead"),
                     control = list(maxit = 2000))
        if(fit$convergence == 10){
            t <- 1
            maxt <- 10
            while(t <= maxt){
                fit <- optim(par = lapply(fit$par, FUN = function(X) rnorm(n=1, mean=X, sd=0.01)),
                             fn = nll,
                             data = as.list(cfdata),
                             disp.mod = paste(m.list[i]),
                             method = c("Nelder-Mead"),
                             control = list(maxit = 2000))
                if(fit$convergence == 0) break
                t <- t + 1
            }
        }
        estimates.temp <- cbind(data.frame(Species = c("CF"),
                                           Model = m.list[i],
                                           Parameters = length(fit$par),
                                           negloglike = fit$value,
                                           set = j,
                                           converge = fit$convergence),
                                t(as.data.frame(fit$par)))
        estimates <- rbind.fill(estimates, estimates.temp)
        startingvalues.temp <- as.list(fit$par)
    }
}

# Summary
estimates.summary <- subset(estimates, Model == c("model0"))
for(i in 2:length(m.list)){
    temp1 <- subset(estimates, Model == m.list[i] & Species == "CS")
    temp2 <- subset(estimates, Model == m.list[i] & Species == "CF")
    estimates.summary <- rbind(estimates.summary, head(temp1[order(temp1$negloglike), ], n = 1))
    estimates.summary <- rbind(estimates.summary, head(temp2[order(temp2$negloglike), ], n = 1))
}
estimates.summary <- estimates.summary[order(estimates.summary$Species), ]

print(estimates.summary)

# Calculate AIC
estimates.summary$AIC <- 2 * estimates.summary$Parameters + (2 * estimates.summary$negloglike)

# Print to file
tab <- xtable(estimates.summary[, -5], digits = 5)
## tex
print(x = tab, type = "latex", file = "results-fits.tex", digits = 5)
## html
print(x = tab, type = "html", file = "results-fits.html", digits = 5)

