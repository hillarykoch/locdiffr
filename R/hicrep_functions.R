# Variance stabilizing transformation
vstran  <- function(d){
    
    x1r <- rank(d[,1], ties.method = "random")
    x2r <- rank(d[,2], ties.method = "random")
    x1.cdf.func <- ecdf(x1r); x2.cdf.func = ecdf(x2r)
    x1.cdf <- x1.cdf.func(x1r)
    x2.cdf <- x2.cdf.func(x2r)
    new_d <- cbind(x1.cdf, x2.cdf)
    
    return(new_d)
}