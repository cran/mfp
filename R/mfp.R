mfp <- function(formula = formula(data), data = parent.frame(), family = gaussian,
     subset, na.action, init, alpha = 0.05, select = 1, verbose = FALSE, x = TRUE, y = TRUE)
{
# Version 1.2.0     26072004
#
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode="function")
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
cox <- (family$family == "Cox")
if(cox) require(survival)
m <- match.call(expand = FALSE)
temp <- c("", "formula", "data", "weights", "subset", "na.action")
m <- m[match(temp, names(m), nomatch = 0)]
special <- c("strata", "fp")
Terms <- if (missing(data))
        terms(formula, special)
else terms(formula, special, data = data)
    m$formula <- Terms
    m$alpha <- m$select <- m$scale <- m$family <- m$verbose <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Y <- model.extract(m, "response") #
    weights <- model.extract(m, "weights")
#
    offset <- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if (tt == 0 & cox)
        rep(0, nrow(Y))
    else if (tt == 0 & !cox)
        rep(0, length(Y))
    else if (tt == 1)
        m[[offset]]
    else {
        ff <- m[[offset[1]]]
        for (i in 2:tt) ff <- ff + m[[offset[i]]]
        ff
    }
    attr(Terms, "intercept") <- 1
    dropx <- NULL
    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- temp$terms
        if (length(temp$vars) == 1)
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel=TRUE)
        strats <- as.numeric(strata.keep)
    }
    if (length(dropx)) # dropx is number of stratification variables
        newTerms <- Terms[-dropx]  # newTerms are Terms w/o stratif. var
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    if (missing(init))
        init <- NULL
#
# Set up lists for fitter
#
    nx <- ncol(X) - 1
    nobs <- nrow(X)
    df.list <- rep(1, nx)
    scale.list <- rep(FALSE, nx)
    alpha.list <- rep(alpha, nx)
    select.list <- rep(select, nx)
#
# Deal with FP terms
#
    fp.pos <- grep("fp", dimnames(X)[[2]])
    fp.mpos <- attributes(Terms)$specials$fp
    fp.xpos <- unlist(attributes(Terms)$specials) - 1
    if(length(fp.pos) > 0) {
        fp.pos <- fp.pos - 1    # without intercept
        fp.data <- m[, fp.mpos, drop=FALSE]
        df.list[fp.pos] <- unlist(lapply(fp.data, attr, "df"))
        scale.list[fp.pos] <- unlist(lapply(fp.data, attr, "scale"))
        alpha.list[fp.pos] <- unlist(lapply(fp.data, attr, "alpha"))
        alpha.list[sapply(alpha.list, is.na)] <- alpha
        select.list[fp.pos] <- unlist(lapply(fp.data, attr, "select"))
        select.list[sapply(select.list, is.na)] <- select
        names <- dimnames(X)[[2]]
        names[fp.pos + 1] <- unlist(lapply(fp.data, attr, "name"))
        xnames <- names[-1]
        xnames[-fp.pos] <- attr(Terms,"term.labels")[-fp.xpos]
        dimnames(X)[[2]] <- names
    }
    unlist(lapply(m, attr, "name"))
#
# Cox fit
#
    if(cox) {
        if(!inherits(Y, "Surv"))
            stop("Response must be a survival object")
        type <- attr(Y, "type")
        if(type != "right")
            stop("The data must be right censored")
        X <- X[, -1, drop=FALSE]  # remove intercept
        control <- coxph.control()
        method <- "efron"
        fit <- mfp.fit(X, Y, TRUE, FALSE, df.list, scale.list, alpha.list, select.list,
                  verbose = verbose, strata=strats, offset=offset, init, control,
                  weights = weights, method = method, rownames = row.names(m),
                  xnames = xnames)
        if(is.character(fit)) {
            fit <- list(fail = fit)
            attr(fit, "class") <- c("mfp", "coxph")
        }
        else attr(fit, "class") <- c("mfp", fit$method)
        fit$n <- nobs
        fit$family <- family
    }
    else {
#
# GLM fit
#
        gauss <- (family$family == "gaussian")
        fit <- mfp.fit(X, Y, FALSE, gauss, df.list, scale.list, alpha.list, select.list,
                  verbose = verbose, family = family, xnames = xnames)
        attr(fit, "class") <- c("mfp", "glm", "lm")
    }
# glm add-on (to derive var)
   if(!cox) {
        dispersion <-
        if (any(fit$family$family == c("poisson","binomial")))
            1
        else if (fit$df.residual > 0) {
            if (any(fit$weights == 0))
                warning("observations with zero weight ", "not used for calculating dispersion")
            sum(fit$weights * fit$residuals^2)/fit$df.residual
        } else Inf
   p <- fit$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- fit$qr
        aliased <- is.na(coef(fit))
        coef.p <- fit$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        fit$var <- dispersion * covmat.unscaled
    }
   }
# back-transformation of coefs and vars
  if(length(fit$df.final)>1) {
   i2 <- if(cox) 0 else 1
   tp <- if(cox) double() else 1
   for(i in seq(fit$df.final)) {
    if(fit$df.final[i]!=0) {
    p <- 1; if(fit$df.final[i]==4) p <- 2
    i1 <- i2+1; i2 <- i2+p
    tp <- c(tp, fit$scale[i,2]^fit$powers[i,1:p])
    fit$coefficients[i1:i2] <- fit$coefficients[i1:i2]/(fit$scale[i,2]^fit$powers[i,1:p])  # coefs
    }
   }
  if(length(tp)>0) fit$var <- fit$var/(tp%*%t(tp))  # vars
  } else {
    if(fit$df.final!=0) {
      p <- 1; if(fit$df.final==4) p <- 2
      if(cox) {
       tp <- as.vector(fit$scale[2]^fit$powers[1:p])
       fit$coefficients <- fit$coefficients/(fit$scale[2]^fit$powers[1:p])  # coefs
       fit$var <- fit$var/(tp%*%t(tp))  # vars
      } else {
       tp <- c(1,as.vector(fit$scale[2]^fit$powers[1:p]))
       fit$coefficients[-1] <- fit$coefficients[-1]/(fit$scale[2]^fit$powers[1:p])  # coefs
       fit$var <- fit$var/(tp%*%t(tp))  # vars
      }
    }
  }
# Cox add-on (needed for summary.coxph)
    if(cox) {
      if(length(fit$coefficients) && is.null(fit$wald.test)) {
        nabeta <- !is.na(fit$coefficients)
        if(is.null(init))
          temp <- fit$coefficients[nabeta]
        else temp <- (fit$coefficients - init)[nabeta]
      }
    if(exists("coxph.wtest")) tester <- get("coxph.wtest")
    else tester <- getFromNamespace("coxph.wtest","survival")

    fit$wald.test <- tester(fit$var[nabeta, nabeta],
                temp, .Machine$double.eps^0.75)$test
    }
#
    if (x) fit$X <- X
    if (y) fit$y <- Y
    fit$terms <- Terms
    fit$call <- call
    fit
}
