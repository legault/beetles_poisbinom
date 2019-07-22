model12 <- function(patch1.before, patch2.before, patch1.other.before, p0, b1, b2, b3, ...){
    lmod <- p0 + b1 * patch1.before + b2 * patch2.before + b3 * patch1.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
