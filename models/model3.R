model3 <- function(patch1.before, patch1.other.before, p0, b1, b3, ...){
    lmod <- p0 + b1 * patch1.before + b3 * patch1.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
