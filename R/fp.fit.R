fp.fit <- function(X, Y, df, dfr, cox, gauss, shift, scale, ...)
{
#
# Version 1.1.1     09 dec 03
#
    X <- X[,apply(X, 2, function(x) !all(x==0)), drop=FALSE] # to avoid numerical problems
    x <- X[, 1] # x is first column of X                 # induced by fp.gen if df=0
    Xcov <- X[, -1, drop=FALSE]   # covariates
    ncov <- ncol(Xcov)
    nobs <- nrow(X)
    pwrs <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
    npwrs <- length(pwrs)
    dev2 <- dev4 <- 1000000000000000    
#
# Set up fitters
#
    if(cox) {
        if(exists("coxph.fit")) fitter <- get("coxph.fit")
        else fitter <- getFromNamespace("coxph.fit","survival")
        dv <- "loglik"
    }
    else {
        fitter <- get("glm.fit")
        dv <- "deviance"
    }
    n <- 1 + cox
#
# Null and linear models
#
    dev <- (-2)^cox * fitter(X, Y, ...)[[dv]]
    dev1 <- dev[n]
    if(!ncov)
        dev0 <- dev[1]
    else dev0 <- (-2)^cox * fitter(Xcov, Y, ...)[[dv]][n]
    if(df > 1) {
        for(i in 1:npwrs) {
#
# Find best single power transformation
#
            x.fp <- fp.gen(x, pwrs[i], shift, scale)
            dev <- (-2)^cox * fitter(cbind(x.fp, Xcov), Y, ...)[[dv
                ]][n]
            if(!is.null(dev) & dev < dev2) {
                dev2 <- dev
                pwr2 <- pwrs[i]
            }
            if(df == 4) {
#
# Find best two power transformation
#
                x.fp <- fp.gen(x, c(pwrs[i], pwrs[i]), shift, 
                  scale)
                dev <- (-2)^cox * fitter(cbind(x.fp, Xcov), Y, 
                  ...)[[dv]][n]
                if(!is.null(dev) & dev < dev4) {
                  dev4 <- dev
                  pwr4 <- c(pwrs[i], pwrs[i])
                }
                j <- i + 1
                while(j <= npwrs) {
                  x.fp <- fp.gen(x, c(pwrs[i], pwrs[j]), shift, scale)
                  dev <- (-2)^cox * fitter(cbind(x.fp, Xcov), Y,
                    ...)[[dv]][n]
                  if(!is.null(dev) & dev < dev4) {
                    dev4 <- dev
                    pwr4 <- c(pwrs[i], pwrs[j])
                  }
                  j <- j + 1
                }
            }
        }
    }
#
# Output
#   
    dev <- c(dev0, dev1, dev2, dev4)
    if(df < 4)
        dev[4] <- pwr4 <- NA
    if(df < 2)
        dev[3] <- pwr2 <- NA
    if(gauss)
        dev <- nobs * (1 + log((2 * pi * dev)/nobs))
    fit <- list(pwr4 = pwr4, pwr2 = pwr2, dev4 = dev[4], dev2 = dev[3], 
        dev1 = dev[2], dev0 = dev[1], nobs = nobs, dfr = dfr, df = df, 
        gauss = gauss)
    return(fit)
}
