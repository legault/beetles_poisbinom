model2 <- function(patch1.before, patch2.before, p0, b1, b2, ...){
    lmod <- p0 + b1 * patch1.before + b2 * patch2.before
    return(exp(lmod) / (exp(lmod) + 1))
}
