model10 <- function(patch2.before, p0, b2, ...){
    lmod <- p0 + b2 * patch2.before
    return(exp(lmod) / (exp(lmod) + 1))
}
