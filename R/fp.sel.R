fp.sel <- function(fit, alpha = 0.05, select = 1)
{
#
# Version 1.3     27.03.2005
#
# Calculate deviance differences & p-values
#
#    dfr <- fit$dfr - fit$df # residual df     # das stimmt doch nicht

# scale: dispersion in fit einfuegen!!!
    dd.null <- fit$dev0 - min(fit$dev1, fit$dev2, fit$dev4, na.rm=TRUE)  # test drop1 indep of FP
    p.null <- 1 - pchisq(dd.null, fit$df)
    if(fit$df > 1) {
        dd.lin <- fit$dev1 - min(fit$dev2, fit$dev4, na.rm=TRUE)
        p.lin <- 1 - pchisq(dd.lin, fit$df - 1)
        if(fit$df > 2) {
            dd.FP <- fit$dev2 - fit$dev4
            p.FP <- 1 - pchisq(dd.FP, 2)
        }
        else p.FP <- NA
    }
    else {
        p.lin <- NA
        p.FP <- NA
    }
    if(p.null > select) {
        df <- 0
        pwrs <- c(NA, NA)
        dev <- fit$dev0
    }
    else {
        if(fit$df > 1) {
            if(p.lin > alpha)
                df <- 1
            else {
                if(fit$df > 2) {
                  if(p.FP > alpha)
                    df <- 2
                  else {
                    df <- 4
                    pwrs <- fit$pwr4
                    dev <- fit$dev4
                  }
                }
                else df <- 2
            }
        }
        else df <- 1
    }
    if(df == 1) {
        pwrs <- c(1, NA)
        dev <- fit$dev1
    }
    if(df == 2) {
        pwrs <- c(fit$pwr2, NA)
        dev <- fit$dev2
    }
    results <- list(p.null = p.null, p.lin = p.lin, p.FP = p.FP, df = df, 
        pwrs = pwrs, dev = dev)
    return(list(results=results, fit=fit))
}
