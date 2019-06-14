model8 <- function(patch2.other.before, p0, b4, ...){
    lmod <- p0 + b4 * patch2.other.before
    return(exp(lmod) / (exp(lmod) + 1))
}
