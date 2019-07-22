model13 <- function(patch2.before, patch2.other.before, p0, b2, b4, ...){
    lmod <- p0 + b2 * patch2.before + b4 * patch2.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
