model5alt <- function(patch1.before, patch2.before, patch1.other.before, patch2.other.before, p0, b1, b2, b3, b4, ...){
    e1 <- (patch1.before + patch2.before) / 2
    delta <- patch1.before - e1
    e2 <- (patch1.other.before + patch2.other.before) / 2
    delta.other <- patch1.other.before - e2
    lmod <- p0 + b1 * e1 + b2 * delta + b3 * e2 + b4 * delta.other
    return(exp(lmod) / (exp(lmod) + 1))
}
