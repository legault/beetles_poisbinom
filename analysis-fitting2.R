# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently
# compute the exact pmf of a Poisson binomial distribution with function 'dpoibin()'
library(poibin)
# Source plyr for rbind.full function
library(plyr)
# Source xtable for saving tables
library(xtable)
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

# Set seed
set.seed(20190422)


# Rounding function
mround <- function(x, base){
        base * round(x / base)
}

# Source data
cs <- read.csv("data/csDat.csv", header = TRUE)
# Round counts
cs$numCS1b.r <- mround(cs$numCS1b, 10)
cs$numCS2b.r <- mround(cs$numCS2b, 10)
cs$numCF1b.r <- mround(cs$numCF1b, 10)
cs$numCF2b.r <- mround(cs$numCF2b, 10)
cs$numCS1a.r <- mround(cs$numCS1a, 10)
cs$numCS2a.r <- mround(cs$numCS2a, 10)
cs$numCF1a.r <- mround(cs$numCF1a, 10)
cs$numCF2a.r <- mround(cs$numCF2a, 10)
# Sum patches
cs$totCSb.r <- cs$numCS1b.r + cs$numCS2b.r
cs$totCSa.r <- cs$numCS1a.r + cs$numCS2a.r
cs$totCFb.r <- cs$numCF1b.r + cs$numCF2b.r
cs$totCFa.r <- cs$numCF1a.r + cs$numCF2a.r
## Totals
cs$totCS.r <- cs$totCSa.r + cs$totCSb.r
cs$totCF.r <- cs$totCFa.r + cs$totCFb.r
## Diff
cs$totCSdiff.r <- cs$totCSa.r - cs$totCSb.r
cs$totCFdiff.r <- cs$totCFa.r - cs$totCFb.r


## ## cs1 <- subset(cs, gen == 1 & totCSdiff.r == 0 & totCS.r > 0)
## ## nrow(cs1)
## ## cs2 <- subset(cs, gen == 2 & totCSdiff.r == 0 & totCS.r > 0)
## ## nrow(cs2)
## ## cs3 <- subset(cs, gen == 3 & totCSdiff.r == 0 & totCS.r > 0)
## ## nrow(cs3)

## cs1 <- subset(cs, gen == 1 & totCSdiff.r == 0 & totCFdiff.r == 0 & totCS.r > 0 & totCF.r > 0)
## nrow(cs1)
## cs2 <- subset(cs, gen == 2 & totCSdiff.r == 0 & totCFdiff.r == 0 & totCS.r > 0 & totCF.r > 0)
## nrow(cs2)
## cs3 <- subset(cs, gen == 3 & totCSdiff.r == 0 & totCFdiff.r == 0 & totCS.r > 0 & totCF.r > 0)
## nrow(cs3)

# Extract relevant information from data
## CS
## Subset
cs.r <- subset(cs, totCSdiff.r == 0 & totCFdiff.r == 0 & totCS.r > 0)
## Create smaller table of relevant data
csdata <- data.frame(patch1.before = cs.r[, "numCS1b.r"],
                     patch2.before = cs.r[, "numCS2b.r"],
                     patch1.after = cs.r[, "numCS1a.r"],
                     patch2.after = cs.r[, "numCS2a.r"],
                     patch1.other.before = cs.r[, "numCF1b.r"],
                     patch2.other.before = cs.r[, "numCF2b.r"])
## CF
cf.r <- subset(cs, totCSdiff.r == 0 & totCFdiff.r == 0 & totCF.r > 0)
## Create smaller table of relevant data
cfdata <- data.frame(patch1.before = cf.r[, "numCF1b.r"],
                     patch2.before = cf.r[, "numCF2b.r"],
                     patch1.after = cf.r[, "numCF1a.r"],
                     patch2.after = cf.r[, "numCF2a.r"],
                     patch1.other.before = cf.r[, "numCS1b.r"],
                     patch2.other.before = cf.r[, "numCS2b.r"])


# Source negative log-likelihood function (general)
source("nll.R")$value
# Create a list models of models to optimize
m.list <- grep("model", ls(), value = TRUE)
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
### Estimate baseline dispersal rate (p0)
cs.gen1 <- subset(cs.r, gen == 1 & numCS1b.r > 0) # generation 1
propd <- mean(cs.gen1$numCS2a.r) / mean(cs.gen1$numCS1b.r) # proportion dispersed
### Set estimate propd as starting value for p0
startingvalues <- list(p0 = log(propd / (1 - propd)))
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
numr <- 10
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
            fit <- optim(par = as.list(fit$par),
                     fn = nll,
                     data = as.list(csdata),
                     disp.mod = paste(m.list[i]),
                     method = c("Nelder-Mead"),
                     control = list(maxit = 2000))
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
### Estimate baseline dispersal rate (p0)
cf.gen1 <- subset(cf.r, gen == 1 & numCF2b.r > 0) # generation 1
propd <- mean(cf.gen1$numCF1a.r) / mean(cf.gen1$numCF2b.r) # proportion dispersed
### Set estimate propd as starting value for p0
startingvalues <- list(p0 = log(propd / (1 - propd)))
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
numr <- 10
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
            fit <- optim(par = as.list(fit$par),
                     fn = nll,
                     data = as.list(cfdata),
                     disp.mod = paste(m.list[i]),
                     method = c("Nelder-Mead"),
                     control = list(maxit = 2000))
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

tab <- xtable(estimates.summary[, -5], digits = 4)
print(x = tab, type = "html", file = "results-fits2.html", digits = 4)

