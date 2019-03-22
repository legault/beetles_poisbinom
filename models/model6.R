model6 <- function(patch1.other.before, p0, b3, ...){
    lmod <- p0 + b3 * patch1.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
