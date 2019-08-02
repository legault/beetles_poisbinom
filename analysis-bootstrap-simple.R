# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently
# compute the exact pmf of a Poisson binomial distribution with function 'dpoibin()'
library(poibin)
# Source plyr for rbind.full function
library(plyr)
# Source xtable for saving tables
library(xtable)
# Source best-fit models
source("models/model6.R")$value # CS and CF

# Set seed
set.seed(20190719)

### Number of bootstraps
numr <- 5000

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

# Extract relevant information from data
## CS
## Subset
cs.r <- subset(cs, totCSdiff.r == 0 & totCFdiff.r == 0 & totCS.r > 0)
## Create smaller table of relevant data
csdata <- data.frame(gen = cs.r[, "gen"],
                     patch1.before = cs.r[, "numCS1b.r"],
                     patch2.before = cs.r[, "numCS2b.r"],
                     patch1.after = cs.r[, "numCS1a.r"],
                     patch2.after = cs.r[, "numCS2a.r"],
                     patch1.other.before = cs.r[, "numCF1b.r"],
                     patch2.other.before = cs.r[, "numCF2b.r"])

## CF
cf.r <- subset(cs, totCSdiff.r == 0 & totCFdiff.r == 0 & totCF.r > 0)
## Create smaller table of relevant data
cfdata <- data.frame(gen = cf.r[, "gen"],
                     patch1.before = cf.r[, "numCF1b.r"],
                     patch2.before = cf.r[, "numCF2b.r"],
                     patch1.after = cf.r[, "numCF1a.r"],
                     patch2.after = cf.r[, "numCF2a.r"],
                     patch1.other.before = cf.r[, "numCS1b.r"],
                     patch2.other.before = cf.r[, "numCS2b.r"])

# Source negative log-likelihood function (general)
source("nll.R")$value
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
### Set estimates as starting values for p0
startingvalues <- list(p0 = -1.10995,
                       b3 = 0.00231)
### Run optimizer
for(j in 1:numr){
    csdata.gen1 <- subset(csdata, gen == 1)
    csdata.gen2 <- subset(csdata, gen == 2)
    csdata.gen3 <- subset(csdata, gen == 3)
    rands.gen1 <- sample(1:nrow(csdata.gen1), size = nrow(csdata.gen1), replace = TRUE)
    rands.gen2 <- sample(1:nrow(csdata.gen2), size = nrow(csdata.gen2), replace = TRUE)
    rands.gen3 <- sample(1:nrow(csdata.gen3), size = nrow(csdata.gen3), replace = TRUE)
    csdata.gen1 <- csdata.gen1[rands.gen1, ]
    csdata.gen2 <- csdata.gen2[rands.gen2, ]
    csdata.gen3 <- csdata.gen3[rands.gen3, ]
    data.temp <- rbind(csdata.gen1, csdata.gen2, csdata.gen3)[, -1]
    fit <- optim(par = startingvalues,
                 fn = nll,
                 data = as.list(data.temp),
                 disp.mod = model6,
                 method = c("Nelder-Mead"),
                 control = list(maxit = 2000))
    if(fit$convergence == 10){
        t <- 1
        maxt <- 10
        while(t <= maxt){
            fit <- optim(par = lapply(fit$par, FUN = function(X) rnorm(n=1, mean=X, sd=0.01)),
                         fn = nll,
                         data = as.list(data.temp),
                         disp.mod = model6,
                         method = c("Nelder-Mead"),
                         control = list(maxit = 2000))
            if(fit$convergence == 0) break
            t <- t + 1
        }
    }
    estimates.temp <- cbind(data.frame(Species = c("CS"),
                                       Model = c("model6"),
                                       Parameters = length(fit$par),
                                       negloglike = fit$value,
                                       set = j,
                                       converge = fit$convergence),
                            t(as.data.frame(fit$par)))
    estimates <- rbind.fill(estimates, estimates.temp)
}

# CF
startingvalues <- list(p0 = -1.93904,
                       b3 = -0.0066)
### Run optimizer
for(j in 1:numr){
    cfdata.gen1 <- subset(cfdata, gen == 1)
    cfdata.gen2 <- subset(cfdata, gen == 2)
    cfdata.gen3 <- subset(cfdata, gen == 3)
    rands.gen1 <- sample(1:nrow(cfdata.gen1), size = nrow(cfdata.gen1), replace = TRUE)
    rands.gen2 <- sample(1:nrow(cfdata.gen2), size = nrow(cfdata.gen2), replace = TRUE)
    rands.gen3 <- sample(1:nrow(cfdata.gen3), size = nrow(cfdata.gen3), replace = TRUE)
    cfdata.gen1 <- cfdata.gen1[rands.gen1, ]
    cfdata.gen2 <- cfdata.gen2[rands.gen2, ]
    cfdata.gen3 <- cfdata.gen3[rands.gen3, ]
    data.temp <- rbind(cfdata.gen1, cfdata.gen2, cfdata.gen3)[, -1]
    fit <- optim(par = startingvalues,
                 fn = nll,
                 data = as.list(data.temp),
                 disp.mod = model6,
                 method = c("Nelder-Mead"),
                 control = list(maxit = 2000))
    if(fit$convergence == 10){
        t <- 1
        maxt <- 10
        while(t <= maxt){
            fit <- optim(par = lapply(fit$par, FUN = function(X) rnorm(n=1, mean=X, sd=0.01)),
                         fn = nll,
                         data = as.list(data.temp),
                         disp.mod = model6,
                         method = c("Nelder-Mead"),
                         control = list(maxit = 2000))
            if(fit$convergence == 0) break
            t <- t + 1
        }
    }
    estimates.temp <- cbind(data.frame(Species = c("CF"),
                                       Model = c("model6"),
                                       Parameters = length(fit$par),
                                       negloglike = fit$value,
                                       set = j,
                                       converge = fit$convergence),
                            t(as.data.frame(fit$par)))
    estimates <- rbind.fill(estimates, estimates.temp)
}

## # Summary
## tab <- xtable(estimates, digits = 5)
## print(x = tab, type = "latex", file = "results-bootstrap.tex", digits = 4)

write.table(estimates, file="results-bootstrap-simple.csv", sep = ",", quote = FALSE, row.names = FALSE)
