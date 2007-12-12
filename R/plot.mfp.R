plot.mfp <- function (x, var=NULL, ref.zero=TRUE, ask=TRUE, ...) 
{
    if (!inherits(x, "mfp")) 
        stop("This is not an mfp object")
    name <- dimnames(x$powers)[[1]]
    choices <- name
    if(is.null(var)) {
      w <- which(is.na(x$powers[,1]))
      if(length(w)==0) { pick <- seq(name) } 
      else { pick <- seq(name)[-w] }
    } else { 
      pick <- which(name %in% var)
    }
    int <- as.numeric(x$family[["family"]] != "Cox")
    for(ip in pick) {
        namex <- name[ip]
        if (is.null(x$X)) 
            stop("you did not specify x=TRUE in the fit")
        if (any(dimnames(x$X)[[2]] == namex, na.rm = TRUE)) {
            tmpx <- x$X[, namex]
            ix <- which(dimnames(x$X)[[2]] == namex)
        }
        else {
            tmpx <- eval(as.name(namex))
        }
        ord <- order(tmpx)
    	tmpx <- sort(tmpx)
#
        npwrsx <- sum(is.finite(x$powers[ip, ]))
        if (npwrsx > 0) {
            if (ip > 1) 
                posx <- int + sum(is.finite(x$powers[seq(ip-1), ])) + seq(npwrsx)
            else posx <- int + seq(npwrsx)
			xtmp <- t(matrix(tmpx,ncol=length(tmpx),nrow=npwrsx,byrow=TRUE)^x$powers[ip,1:npwrsx])
            fx <- xtmp %*% x$coef[posx]
            if (ref.zero) {
             xcenter <- mean(tmpx, na.rm=TRUE)^x$powers[ip,1:npwrsx]
             ycenter <- as.vector(xcenter %*% x$coef[posx])
             fx <- fx - ycenter
            }
        }

# Plots
   if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
#
        if (int) {                      # generalized linear model
            if (npwrsx > 0) {
                plot(tmpx, fx, xlab = namex, ylab = paste("Linear predictor", 
                  sep = ""), type = "l", ...)
                pres <- x$residuals[ord] + fx
                plot(tmpx, pres, xlab = namex, ylab = "Partial residuals", 
                  ...)
				ok <- is.finite(tmpx) & is.finite(pres)
                fl <- lowess(tmpx[ok], pres[ok], iter = 0)
                lines(fl$x, fl$y, lwd = 1, col = "red")
            }
        }
        else {                          # Cox proportional hazards model
            require(survival)
# Martingale residuals
            x0 <- coxph(x$y ~ 1)
            res0 <- resid(x0, type = "mart")[ord]
            plot(tmpx, res0, xlab = namex, ylab = "Martingale residuals", 
                type = "p", ...)
			ok <- is.finite(tmpx) & is.finite(res0)
			fl <- lowess(tmpx[ok], res0[ok], iter = 0)
            lines(fl$x, fl$y, lwd = 1, col = "red")
# Fitted function
            if (npwrsx > 0) {
                plot(tmpx, fx, xlab = namex, ylab = "Linear predictor", 
                  type = "l", ...)
                pres <- x$residuals[ord] + fx
                plot(tmpx, pres, xlab = namex, ylab = "Partial residuals", 
                  ...)
				ok <- is.finite(tmpx) & is.finite(pres)
                fl <- lowess(tmpx[ok], pres[ok], iter = 0)
                lines(fl$x, fl$y, lwd = 1, col = "red")
            }
        }
    }
}
