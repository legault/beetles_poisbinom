# Source 'poibin' package, which uses the Discrete Fourier Transform (DFT) to efficiently
# compute the exact pmf of a Poisson binomial distribution with function 'dpoibin()'
library(poibin)
# Source plyr for rbind.full function
library(plyr)
# Source xtable for saving tables
library(xtable)

# Source models
source("models/model1.R")$value
source("models/model2.R")$value
source("models/model3.R")$value
source("models/model4.R")$value
source("models/model5.R")$value
source("models/model6.R")$value
source("models/model7.R")$value
# Source data
cs <- read.csv("data/csDat.csv", header = TRUE)
# Source negative log-likelihood function (general)
source("nll.R")$value

# Extract relevant information from data
cs$totCSb <- cs$numCS1b + cs$numCS2b # total CS in both patches (before)
cs$totCSa <- cs$numCS1a + cs$numCS2a # total CS in both patches (after)
cs$totCSdiff <- cs$totCSb - cs$totCSa # calculate difference between before and after
cs1 <- subset(cs, gen == 1 & totCSdiff == 0 & totCSb > 0) # keep only data for gen1 and diff 0
## Create smaller table of relevant data
cs1data <- data.frame(patch1.before = cs1[, "numCS1b"],
                      patch2.before = cs1[, "numCS2b"],
                      patch1.after = cs1[, "numCS1a"],
                      patch2.after = cs1[, "numCS2a"],
                      patch1.other.before = cs1[, "numCF1b"],
                      patch2.other.before = cs1[, "numCF2b"])

# Calculate proportion dispered
propd <- mean(cs1data$patch2.after / cs1data$patch1.before)
# Set starting values, with propd baseline estimate for p0; assume other parameters have small effect
set.seed(20190422)
startingvalues <- data.frame(p0 = log(propd / (1 - propd)),
                             b1 = rnorm(n = 1, mean = 0, sd = .1),
                             b2 = rnorm(n = 1, mean = 0, sd = .1),
                             b3 = rnorm(n = 1, mean = 0, sd = .1),
                             b4 = rnorm(n = 1, mean = 0, sd = .1))

# Create a list models of models to optimize
m.list <- grep("model", ls(), value = TRUE)
# Create empty table of estimates
estimates <- data.frame(Model = c(),
                        Parameters = c(),
                        negloglike = c(),
                        p0 = c(),
                        b1 = c(),
                        b2 = c(),
                        b3 = c(),
                        b4 = c(),
                        converge = c())
# Optimize
## Run optimizer
for(i in 1:7){
    arg.temp <- formalArgs(paste(m.list[i]))
    startingvalues.temp <- as.list(startingvalues[, which(names(startingvalues) %in% arg.temp)])
    fit <- optim(par = startingvalues.temp, # use true values as initial values
                 fn = nll,
                 data = as.list(cs1data),
                 disp.mod = paste(m.list[i]),
                 control = list(maxit = 1000))
    estimates.temp <- cbind(data.frame(Model = m.list[i],
                                       Parameters = length(fit$par),
                                       negloglike = fit$value,
                                       converge = fit$convergence),
                            t(as.data.frame(fit$par)))
    estimates <- rbind.fill(estimates, estimates.temp)
}

# Save results
tab <- xtable(estimates, digits = 4)
print(x = tab, type = "html", file = "results-fits.html", digits = 4)

