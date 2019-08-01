
# Source xtable for saving tables
library(xtable)

results <- read.csv(file = "results-bootstrap-all.csv", header = TRUE)

cs <- subset(results, Species == "CS" & converge == 0)
cs.summary <- data.frame(Species = "CS",
                         Model = "model14",
                         mean.p0 = mean(cs$p0, na.rm = TRUE),
                         mean.b1 = mean(cs$b1, na.rm = TRUE),
                         mean.b2 = mean(cs$b2, na.rm = TRUE),
                         mean.b4 = mean(cs$b4, na.rm = TRUE),
                         lower95.p0 = quantile(cs$p0, probs = c(0.025)),
                         lower95.b1 = quantile(cs$b1, probs = c(0.025)),
                         lower95.b2 = quantile(cs$b2, probs = c(0.025)),
                         lower95.b4 = quantile(cs$b4, probs = c(0.025)),
                         upper95.p0 = quantile(cs$p0, probs = c(0.975)),
                         upper95.b1 = quantile(cs$b1, probs = c(0.975)),
                         upper95.b2 = quantile(cs$b2, probs = c(0.975)),
                         upper95.b4 = quantile(cs$b4, probs = c(0.975)),
                         row.names = NULL)

cf <- subset(results, Species == "CF" & converge == 0)
cf.summary <- data.frame(Species = "CF",
                         Model = "model11",
                         mean.p0 = mean(cf$p0, na.rm = TRUE),
                         mean.b2 = mean(cf$b2, na.rm = TRUE),
                         mean.b3 = mean(cf$b3, na.rm = TRUE),
                         lower95.p0 = quantile(cf$p0, probs = c(0.025)),
                         lower95.b2 = quantile(cf$b2, probs = c(0.025)),
                         lower95.b3 = quantile(cf$b3, probs = c(0.025)),
                         upper95.p0 = quantile(cf$p0, probs = c(0.975)),
                         upper95.b2 = quantile(cf$b2, probs = c(0.975)),
                         upper95.b3 = quantile(cf$b3, probs = c(0.975)),
                         row.names = NULL)

write.table(cs.summary, file = "bootCS.csv", sep = ",", row.names = FALSE, quote = FALSE)

write.table(cf.summary, file = "bootCF.csv", sep = ",", row.names = FALSE, quote = FALSE)



                         
                         
