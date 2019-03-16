nll <- function(param, data, disp.mod){
    d1 <- do.call(disp.mod, append(data, param)) # probability disperse from patch1
    d2 <- do.call(disp.mod, append(data, param)) # probability disperse from patch2
    loglike <- dim(length(data$patch1.before))
    for(i in 1:length(data$patch1.before)){
        loglike[i] <- log(dpoibin(kk = data$patch1.after[i],
                                  pp = c(rep(1 - d1[i], data$patch1.before[i]),
                                         rep(d2[i], data$patch2.before[i]))))
    }
    loglike[is.na(loglike)] <- -10000 # punish Inf
    return(-sum(loglike))
}
