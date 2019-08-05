# Create rounding function
mround <- function(x, base){
        base * round(x / base)
}

# Source data
cs <- read.csv("../data/csDat.csv", header = TRUE)
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
csdata <- data.frame(patch1.before = cs.r[, "numCS1b"],
                     patch2.before = cs.r[, "numCS2b"],
                     patch1.after = cs.r[, "numCS1a"],
                     patch2.after = cs.r[, "numCS2a"],
                     patch1.other.before = cs.r[, "numCF1b"],
                     patch2.other.before = cs.r[, "numCF2b"])
## CF
cf.r <- subset(cs, totCSdiff.r == 0 & totCFdiff.r == 0 & totCF.r > 0)
## Create smaller table of relevant data
cfdata <- data.frame(patch1.before = cf.r[, "numCF1b"],
                     patch2.before = cf.r[, "numCF2b"],
                     patch1.after = cf.r[, "numCF1a"],
                     patch2.after = cf.r[, "numCF2a"],
                     patch1.other.before = cf.r[, "numCS1b"],
                     patch2.other.before = cf.r[, "numCS2b"])

# Source parameters
source("../parameters/paramest-cast.R")$value # T cast
source("../parameters/paramest-conf.R")$value # T conf
# Source models
source("../models/model14.R")$value # T cast
source("../models/model11.R")$value # T conf
# Set number of simulations
sims <- 10000

# CS
## Best fit model
cs.stor <- data.frame(Rep = c(),
                      Mean = c(),
                      Lower95 = c(),
                      Upper95 = c(),
                      Patch2Before = c(),
                      Patch2Other = c())
temp.stor <- dim(sims)
for(i in 1:nrow(csdata)){
    temp.stor <- dim(sims)
    temp <- csdata[i, ]
        for(j in 1:sims){
        d1 <- rbinom(n = 1, size = temp[, 1],
                     prob = model14(temp[, 1], temp[, 2], temp[, 6], paramest.cast$p0, paramest.cast$b1, paramest.cast$b2, paramest.cast$b4)) # dispersal patch1
        d2 <- rbinom(n = 1, size = temp[, 2],
                     prob = model14(temp[, 2], temp[, 1], temp[, 5], paramest.cast$p0, paramest.cast$b1, paramest.cast$b2, paramest.cast$b4)) # dispersal patch2
        pred.value <- temp[, 1] - d1 + d2 # abundance after dispersal in patch1
        temp.stor[j] <- temp[, 3] - pred.value
    }
    cs.stor <- rbind(cs.stor, data.frame(Rep = i,
                                         Mean = mean(temp.stor),
                                         Lower95 = quantile(temp.stor, probs = c(0.025))[[1]],
                                         Upper95 = quantile(temp.stor, probs = c(0.975))[[1]],
                                         Patch2Before = temp[, 2],
                                         Patch2Other = temp[, 6]))
}
## Null model
cs.stor0 <- data.frame(Rep = c(),
                      Mean = c(),
                      Lower95 = c(),
                      Upper95 = c())
temp.stor <- dim(sims)
for(i in 1:nrow(csdata)){
    temp.stor <- dim(sims)
    temp <- csdata[i, ]
    for(j in 1:sims){
        d1 <- rbinom(n = 1, size = temp[, 1],
                     prob = model14(temp[, 1], temp[, 2], temp[, 6], -1.07206, 0, 0, 0)) # dispersal patch1
        d2 <- rbinom(n = 1, size = temp[, 2],
                     prob = model14(temp[, 2], temp[, 1], temp[, 5], -1.07206, 0, 0, 0)) # dispersal patch2
        pred.value <- temp[, 1] - d1 + d2 # abundance after dispersal in patch1
        temp.stor[j] <- temp[, 3] - pred.value
    }
    cs.stor0 <- rbind(cs.stor0, data.frame(Rep = i,
                                           Mean = mean(temp.stor),
                                           Lower95 = quantile(temp.stor, probs = c(0.025))[[1]],
                                           Upper95 = quantile(temp.stor, probs = c(0.975))[[1]]))
}
# Divide into relevant sets
setcast1 <- which(cs.stor$Patch2Before > 0 & cs.stor$Patch2Other == 0) # b4 not relevant
setcast2 <- which(cs.stor$Patch2Before > 0 & cs.stor$Patch2Other > 0) # all coefficients relevant
# Reorder according to setcasts
cs.stor <- rbind(cs.stor[setcast1, ], cs.stor[setcast2, ])
cs.stor$Rep <- 1:nrow(cs.stor)
cs.stor0 <- rbind(cs.stor0[setcast1, ], cs.stor0[setcast2, ])
cs.stor0$Rep <- 1:nrow(cs.stor)

# CF
## Best fit model
cf.stor <- data.frame(Rep = c(),
                      Mean = c(),
                      Lower95 = c(),
                      Upper95 = c(),
                      Patch2Before = c(),
                      Patch1Other = c())
temp.stor <- dim(sims)
for(i in 1:nrow(cfdata)){
    temp.stor <- dim(sims)
    temp <- cfdata[i, ]
        for(j in 1:sims){
        d1 <- rbinom(n = 1, size = temp[, 1],
                     prob = model11(temp[, 2], temp[, 5], paramest.conf$p0, paramest.conf$b2, paramest.conf$b3)) # dispersal patch1
        d2 <- rbinom(n = 1, size = temp[, 2],
                     prob = model11(temp[, 1], temp[, 6], paramest.conf$p0, paramest.conf$b2, paramest.conf$b3)) # dispersal patch2
        pred.value <- temp[, 1] - d1 + d2 # abundance after dispersal in patch1
        temp.stor[j] <- temp[, 3] - pred.value
    }
    cf.stor <- rbind(cf.stor, data.frame(Rep = i,
                                         Mean = mean(temp.stor),
                                         Lower95 = quantile(temp.stor, probs = c(0.025))[[1]],
                                         Upper95 = quantile(temp.stor, probs = c(0.975))[[1]],
                                         Patch2Before = temp[, 2],
                                         Patch1Other = temp[, 5]))
}
## Null model
cf.stor0 <- data.frame(Rep = c(),
                      Mean = c(),
                      Lower95 = c(),
                      Upper95 = c())
temp.stor <- dim(sims)
for(i in 1:nrow(cfdata)){
    temp.stor <- dim(sims)
    for(j in 1:sims){
        temp <- cfdata[i, ]
        d1 <- rbinom(n = 1, size = temp[, 1],
                     prob = model11(temp[, 2], temp[, 5], -2.05377, 0, 0)) # dispersal patch1
        d2 <- rbinom(n = 1, size = temp[, 2],
                     prob = model11(temp[, 1], temp[, 6], -2.05377, 0, 0)) # dispersal patch2
        pred.value <- temp[, 1] - d1 + d2 # abundance after dispersal in patch1
        temp.stor[j] <- temp[, 3] - pred.value
    }
    cf.stor0 <- rbind(cf.stor0, data.frame(Rep = i,
                                           Mean = mean(temp.stor),
                                           Lower95 = quantile(temp.stor, probs = c(0.025))[[1]],
                                           Upper95 = quantile(temp.stor, probs = c(0.975))[[1]]))
}

# Divide into relevant sets
setconf1 <- which(cf.stor$Patch2Before == 0 & cf.stor$Patch1Other == 0) # b2 and b3 not relevant
setconf2 <- which(cf.stor$Patch2Before > 0 & cf.stor$Patch1Other == 0) # b3 not relevant
setconf3 <- which(cf.stor$Patch2Before > 0 & cf.stor$Patch1Other > 0) # all coefficients relevant
setconf4 <- which(cf.stor$Patch2Before == 0 & cf.stor$Patch1Other > 0) # b2 not relevant
# Reorder according to setconfs
cf.stor <- rbind(cf.stor[setconf1, ], cf.stor[setconf2, ], cf.stor[setconf3, ], cf.stor[setconf4, ])
cf.stor$Rep <- 1:nrow(cf.stor)
cf.stor0 <- rbind(cf.stor0[setconf1, ], cf.stor0[setconf2, ], cf.stor0[setconf3, ], cf.stor0[setconf4, ])
cf.stor0$Rep <- 1:nrow(cf.stor)

pdf(file = "Goodness.pdf", height = 6, width = 14)
# Plot
transcol1 <- rgb(0, 0, 0, 0.5)
transcol2 <- rgb(0, 0, 1, 0.8)
transcol3 <- rgb(1, 0, 1, 0.8)
transcol4 <- rgb(1, 0, 0, 0.8)
par(mar = c(5, 5, 2, 0), lend = 2, mfrow = c(1, 2), cex.lab = 1.2)
# Panel A
# NUll
plot(NA, xlim = c(1, 40), ylim = c(-60, 60), xlab = "Replicate", ylab = "Observed - Expected", bty = "n", xaxt = "n", yaxt = "n")
axis(1, at = seq(1, 5, 4), tcl = -.5)
axis(1, at = seq(5, 10, 5), tcl = 0)
axis(1, at = seq(10, 40, 5), tcl = -.5)
axis(2, at = seq(-60, 60, 20), tcl = -.5, las = 1)
## Set 1
arrows(x0 = cs.stor0$Rep[1:length(setcast1)], y0 = cs.stor0$Lower95[1:length(setcast1)], y1 = cs.stor0$Upper95[1:length(setcast1)], col = transcol2, angle = 90, length = 0.05, code = 3)
points(x = cs.stor0$Rep[1:length(setcast1)], y = cs.stor0$Mean[1:length(setcast1)], col = transcol2, pch = 15)
## Set 2
shiftx1 <- length(setcast1) + 1
shiftx2 <- length(setcast1) + length(setcast2)
arrows(x0 = cs.stor0$Rep[shiftx1:shiftx2], y0 = cs.stor0$Lower95[shiftx1:shiftx2], y1 = cs.stor0$Upper95[shiftx1:shiftx2], col = transcol3, angle = 90, length = 0.05, code = 3)
points(x = cs.stor0$Rep[shiftx1:shiftx2], y = cs.stor0$Mean[shiftx1:shiftx2], col = transcol3, pch = 15)
## Non-crossing intervals
red1 <- which(cs.stor0$Upper95 < 0 & cs.stor0$Lower95 < 0)
red2 <- which(cs.stor0$Lower95 > 0 & cs.stor0$Upper95 > 0)
redset <- c(red1, red2)
arrows(x0 = redset, y0 = -55, y1 = -60, angle = 45, length = .05, code = 2, col = "black")
## Label
mtext(side = 3, text = expression(paste("(a) Density-independent dispersal (", italic("T. castaneum"), ")", sep = "")), adj = 0, cex = 1.3)
## Legend
legend("topleft", c(expression(paste(italic("N")[1], ", ", italic("N")[2], sep = "")),
                     expression(paste(italic("N")[1], ", ", italic("N")[2], ", ", italic("M")[2], sep = ""))),
       pch = 15, col = c(transcol2, transcol3), inset = 0.05, bty = "n")
# Dotted line
clip(x1 = 1, x2 = 40, y1=-1, y2 = 1)
abline(h = 0, lty = "dotted")
# Panel B
plot(NA, xlim = c(1, 40), ylim = c(-60, 60), xlab = "Replicate", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
axis(1, at = seq(1, 5, 4), tcl = -.5)
axis(1, at = seq(5, 10, 5), tcl = 0)
axis(1, at = seq(10, 40, 5), tcl = -.5)
axis(2, at = seq(-60, 60, 20), tcl = -.5, las = 1)
## Set 1
arrows(x0 = cs.stor$Rep[1:length(setcast1)], y0 = cs.stor$Lower95[1:length(setcast1)], y1 = cs.stor$Upper95[1:length(setcast1)], col = transcol2, angle = 90, length = 0.05, code = 3)
points(x = cs.stor$Rep[1:length(setcast1)], y = cs.stor$Mean[1:length(setcast1)], col = transcol2, pch = 15)
## Set 2
shiftx1 <- length(setcast1) + 1
shiftx2 <- length(setcast1) + length(setcast2)
arrows(x0 = cs.stor$Rep[shiftx1:shiftx2], y0 = cs.stor$Lower95[shiftx1:shiftx2], y1 = cs.stor$Upper95[shiftx1:shiftx2], col = transcol3, angle = 90, length = 0.05, code = 3)
points(x = cs.stor$Rep[shiftx1:shiftx2], y = cs.stor$Mean[shiftx1:shiftx2], col = transcol3, pch = 15)
## Non-crossing intervals
red1 <- which(cs.stor$Upper95 < 0 & cs.stor$Lower95 < 0)
red2 <- which(cs.stor$Lower95 > 0 & cs.stor$Upper95 > 0)
redset <- c(red1, red2)
arrows(x0 = redset, y0 = -55, y1 = -60, angle = 45, length = .05, code = 2, col = "black")
## Label
mtext(side = 3, text = expression(paste("(b) Density-dependent dispersal (", italic("T. castaneum"), ")", sep = "")), adj = 0, cex = 1.3)
# Dotted line
clip(x1 = 1, x2 = 40, y1=-1, y2 = 1)
abline(h = 0, lty = "dotted")
dev.off()
