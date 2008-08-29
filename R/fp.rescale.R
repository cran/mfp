fp.rescale <- function(x)
{
#
cox <- x$family$family == "Cox"
# back-transformation
    x$scaled.coefficients <- x$coefficients
    if (length(x$df.final) > 1) {
        i2 <- if (cox) 0 else 1
        tp <- if (cox) double() else 1
        for(i in seq(x$df.final)) {
            if(x$df.final[i] != 0) {
                p <- 1
                if(x$df.final[i] == 4) p <- 2
                i1 <- i2 + 1
                i2 <- i2 + p
				xs <- x$scale[i, 2]^x$powers[i, 1:p]
                tp <- c(tp, xs)
                x$scaled.coefficients[i1:i2] <- x$coefficients[i1:i2]/xs
            }
        }
        if (length(tp) > 0) 
            x$scaled.var <- x$var/(tp %*% t(tp))
    }
    else {
        if(x$df.final != 0) {
            p <- 1
            if(x$df.final == 4) 
                p <- 2
            if(cox) {
				xs <- x$scale[2]^x$powers[1:p]
                tp <- as.vector(xs)
                x$scaled.coefficients <- x$coefficients/xs
                x$scaled.var <- x$var/(tp %*% t(tp))
            } else {
				xs <- x$scale[2]^x$powers[1:p]
                tp <- c(1, as.vector(xs))
                x$scaled.coefficients[-1] <- x$coefficients[-1]/xs
                x$scaled.var <- x$var/(tp %*% t(tp))
            }
        }
    }
    return(list(coefficients=x$scaled.coefficients, var=x$scaled.var))
}
