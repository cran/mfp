plot.mfp <- function (x, var=NULL, ...) 
{
    if (!inherits(x, "mfp")) 
        stop("This is not an mfp object")
    name <- dimnames(x$powers)[[1]]
    choices <- name
    if(is.null(var)) {
      pick <- seq(name)[-which(is.na(x$powers[,1]))]
    } else { pick <- which(name %in% var) }
    int <- as.numeric(x$family[["family"]] != "Cox")
    for(ip in pick) {
        npwrsx <- sum(!is.na(x$powers[ip, ]))
        if (npwrsx > 0) {
            if (ip > 1) 
                posx <- int + sum(!is.na(x$powers[seq(ip-1), ])) + seq(npwrsx)
            else posx <- int + seq(npwrsx)
            fx <- x$x[, posx, drop = FALSE] %*% x$coef[posx]
        }
        namex <- name[ip]
        if (is.null(x$X)) 
            stop("you did not specify x=T in the fit")
        if (any(dimnames(x$X)[[2]] == namex, na.rm = TRUE)) {
            tmpx <- x$X[, namex]
            ix <- which(dimnames(x$X)[[2]] == namex)
        }
        else {
            tmpx <- eval(as.name(namex))
        }
        ord <- order(tmpx)
        if (int) {
            if (npwrsx > 0) {
                plot(tmpx[ord], fx[ord], xlab = namex, ylab = paste("Linear predictor", 
                  sep = ""), type = "l", ...)
                pres <- x$residuals + fx
                plot(tmpx, pres, xlab = namex, ylab = "Partial residuals", 
                  ...)
                fl <- lowess(tmpx[ord], pres[ord], iter = 0)
                lines(fl$x, fl$y, lwd = 1, col = "red")
            }
        }
        else {
            require(survival)
            x0 <- coxph(x$y ~ 1)
            res0 <- resid(x0, type = "mart")
            plot(tmpx, res0, xlab = namex, ylab = "Martingale residuals", 
                type = "p", ...)
            fl <- lowess(tmpx[ord], res0[ord], iter = 0)
            lines(fl$x, fl$y, lwd = 1, col = "red")
            if (npwrsx > 0) {
                plot(tmpx[ord], fx[ord], xlab = namex, ylab = "Linear predictor", 
                  type = "l", ...)
                pres <- x$residuals + fx
                plot(tmpx, pres, xlab = namex, ylab = "Partial residuals", 
                  ...)
                fl <- lowess(tmpx[ord], pres[ord], iter = 0)
                lines(fl$x, fl$y, lwd = 1, col = "red")
            }
        }
    }
}
