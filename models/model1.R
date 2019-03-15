model1 <- function(patch1.before, p0, b1, ...){
    lmod <- p0 + b1 * patch1.before
    return(exp(lmod) / (exp(lmod) + 1))
}
