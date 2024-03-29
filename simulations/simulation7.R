# Source parameters
source("param.R")$value
# Source model1
source("../models/model7.R")$value
# Set number of simulations
sims <- 100000

# Simulate dispersal
## Set seed
set.seed(20190314)
## Create storage matrix
sim.stor <- matrix(NA, ncol = 6, nrow = sims)
## Fill matrix
sim.stor[, 1] <- rpois(n = sims, lambda = 80) # random abundance before dispersal in patch1
sim.stor[, 2] <- rpois(n = sims, lambda = 50) # random abundance before dispersal in patch2
sim.stor[, 5] <- rpois(n = sims, lambda = 70) # random abundance (other)  before dispersal in patch1
sim.stor[, 6] <- rpois(n = sims, lambda = 20) # random abundance (other) before dispersal in patch2
d1 <- rbinom(n = sims, size = sim.stor[, 1],
             prob = model7(sim.stor[, 5], sim.stor[, 6], param$p0, param$b3, param$b4)) # dispersal patch1
d2 <- rbinom(n = sims, size = sim.stor[, 2],
             prob = model7(sim.stor[, 6], sim.stor[, 5], param$p0, param$b3, param$b4)) # dispersal patch2
sim.stor[, 3] <- sim.stor[, 1] - d1 + d2 # abundance after dispersal in patch1
sim.stor[, 4] <- sim.stor[, 2] - d2 + d1 # abundance after dispersal in patch2
## Rename columns
colnames(sim.stor) <- c("patch1.before", "patch2.before", "patch1.after", "patch2.after", "patch1.other.before", "patch2.other.before")
## Check sims
head(sim.stor)

# Save results
write.table(sim.stor, file = "../data/simulation7.csv", sep = ",", col.names = TRUE, row.names = FALSE)
