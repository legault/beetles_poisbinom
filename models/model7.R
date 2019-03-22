model7 <- function(patch1.other.before, patch2.other.before, p0, b3, b4, ...){
    lmod <- p0 + b3 * patch1.other.before + b4 * patch2.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
