## # Source xtable for saving tables
## library(xtable)

# Create list of models
m.list <- list.files(path = "models", pattern = "*.R")
# Source pmf functions
sapply(m.list, FUN = function(X) source(paste("models/", X, sep = "")))
# Trim ".R" extension from list
m.list <- gsub(pattern = "\\.R$", "", m.list)

paramlist <- c("p0", "b1", "b2", "b3", "b4")

modeltable1 <- t(sapply(1:length(m.list), FUN = function(i){
    paramlist %in% formalArgs(paste(m.list[i]))}))
modeltable1

modeltable2 <- cbind(data.frame(Model = m.list), modeltable1)
colnames(modeltable2)[2:6] <- paramlist
modeltable2
# reorder
modeltable2 <- modeltable2[c(1, 2, 9:16, 3:8), ]
modeltable2

write.table(modeltable2, file = "modeltable.csv", row.names = FALSE, sep = ",")
