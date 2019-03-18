model4 <- function(patch1.before, patch1.other.before, patch2.other.before, p0, b1, b3, b4, ...){
    lmod <- p0 + b1 * patch1.before + b3 * patch1.other.before + b4 * patch2.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
