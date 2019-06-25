nll <- function(param, data, disp.mod){
    # Set probability disperse from patch1
    d1 <- do.call(disp.mod, append(data, param)) # probability disperse from patch1
    if(length(d1) == 1) d1 <- rep(d1, length(data$patch1.before))
    # Switch patch1 and patch2 data
    data.switch <- data # copy data
    data.switch[grep("patch1", names(data), value = TRUE)] <- data[grep("patch2", names(data), value = TRUE)] # switch patch1 with patch2
    data.switch[grep("patch2", names(data), value = TRUE)] <- data[grep("patch1", names(data), value = TRUE)] # switch patch2 with patch1
    # Set probability disperse from patch 2
    d2 <- do.call(disp.mod, append(data.switch, param))
    if(length(d2) == 1) d2 <- rep(d2, length(data$patch1.before))
    loglike <- dim(length(data$patch1.before))
    for(i in 1:length(data$patch1.before)){
        loglike[i] <- log(dpoibin(kk = data$patch1.after[i],
                                  pp = c(1 - d1[i], d2[i]),
                                  wts = c(data$patch1.before[i], data$patch2.before[i])))
    }
    loglike[is.infinite(loglike)] <- -10000 # punish Inf
    return(-sum(loglike))
}
