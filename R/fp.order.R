fp.order <- function(x, y, cox, gauss, xnames, ...)
{
#
# Version 1.1.0     30.12.03
#
# Returns ordering of xvars (not intercept!)
# by Wald test ranking (Cox: by LR test)
#
    int <- as.numeric(!cox)
    nx <- ncol(x)
    if(cox) {
        if(exists("coxph.fit")) fitter <- get("coxph.fit")
        else fitter <- getFromNamespace("coxph.fit","survival")
        fit <- fitter(x, y, ...)
        se <- sqrt(diag(fit$var))
        dev <- -2 * fit$loglik
        p.value <- vector(length=length(xnames))
        for(i in unique(xnames)) {
            ld <- sum(xnames==i)
            ll <- fitter(x[,xnames!=i, drop=FALSE], y, ...)$loglik[2]
            p.value[xnames==i] <- pchisq(2*(fit$loglik[2]-ll), df=ld)
        }
        x.order <- order(p.value[(1 + int):nx])
    }
    else {
        fit <- glm.fit(x, y, ...)
        Iinv <- solve(t(fit$R) %*% fit$R)
        se <- sqrt(fit$dev/fit$df.residual * diag(Iinv))
        dev <- c(fit$null.deviance, fit$deviance)
t.value <- abs(fit$coef/se)
x.order <- rev(order(t.value[(1 + int):nx]))
    }
    if(gauss) {
        nobs <- nrow(x)
        dev <- nobs * (1 + log((2 * pi * dev)/nobs))
    }
    return(list(order = x.order, dev = dev))
}
