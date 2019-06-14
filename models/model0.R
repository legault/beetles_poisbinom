model0 <- function(p0, ...){
    lmod <- p0
    return(exp(lmod) / (exp(lmod) + 1))
}
