mfp <- function (formula = formula(data), data = parent.frame(), family = gaussian, 
    subset = NULL, na.action = na.omit, init = NULL, alpha = 0.05, select = 1, verbose = FALSE, 
    x = TRUE, y = TRUE) 
{
#
# Version 1.3     27.03.2005
#
    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function")
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }
    cox <- (family$family == "Cox")
    if (cox) 
        require(survival)
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
    Y <- model.extract(m, "response")
    weights <- model.extract(m, "weights")
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
if(cox){
    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- temp$terms
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }
}
    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    if (missing(init)) 
        init <- NULL
    nx <- ncol(X) - 1
    nobs <- nrow(X)
    df.list <- rep(1, nx)
    scale.list <- rep(FALSE, nx)
    alpha.list <- rep(alpha, nx)
    select.list <- rep(select, nx)
    fp.pos <- grep("fp", dimnames(X)[[2]])
    fp.mpos <- attributes(Terms)$specials$fp
    if(cox) 
	 fp.xpos <- unlist(attributes(Terms)$specials) - 1
	else
     fp.xpos <- unlist(attributes(Terms)$specials$fp) - 1
	assign <- attr(X, "assign")[-1]
    if (length(fp.pos) > 0) {
        fp.pos <- fp.pos - 1
        fp.data <- m[, fp.mpos, drop = FALSE]
        df.list[fp.pos] <- unlist(lapply(fp.data, attr, "df"))
        scale.list[fp.pos] <- unlist(lapply(fp.data, attr, "scale"))
        alpha.list[fp.pos] <- unlist(lapply(fp.data, attr, "alpha"))
        alpha.list[sapply(alpha.list, is.na)] <- alpha
        select.list[fp.pos] <- unlist(lapply(fp.data, attr, "select"))
        select.list[sapply(select.list, is.na)] <- select
        names <- dimnames(X)[[2]]
        names[fp.pos + 1] <- unlist(lapply(fp.data, attr, "name"))
        xnames <- names[-1]
		tab <- table(assign[-fp.pos])
        xnames[-fp.pos] <- rep(attr(Terms, "term.labels")[-fp.xpos],tab)
        dimnames(X)[[2]] <- names
    } else {
        names <- dimnames(X)[[2]]
        tab <- table(assign)
        xnames <- rep(attr(Terms, "term.labels"),tab)
	}
#    unlist(lapply(m, attr, "name"))
#
# fit cox-mfp model:
#
    if (cox) {
        if (!inherits(Y, "Surv")) 
            stop("Response must be a survival object")
        type <- attr(Y, "type")
        if (type != "right") 
            stop("The data must be right censored")
        X <- X[, -1, drop = FALSE]
        control <- coxph.control()
        method <- "efron"
        fit <- mfp.fit(X, Y, TRUE, FALSE, df.list, scale.list, 
            alpha.list, select.list, verbose = verbose, xnames = xnames,
			strata = strats, offset = offset, init, control, weights = weights, 
            method = method, rownames = row.names(m))
        if (is.character(fit)) {
            fit <- list(fail = fit)
            attr(fit, "class") <- c("mfp", "coxph")
        }
        else attr(fit, "class") <- c("mfp", fit$method)
		attr(fit, "class") <- c("mfp", "coxph")
        fit$n <- nobs
        fit$family <- family
    }
#
# fit glm-mfp model:
#
    else {
        gauss <- (family$family == "gaussian")
        fit <- mfp.fit(X, Y, FALSE, gauss, df.list, scale.list, 
            alpha.list, select.list, verbose = verbose, family = family, 
            xnames = xnames)
        attr(fit, "class") <- c("mfp", "glm", "lm")
    }
# compute dispersion
    if (!cox) {
        dispersion <- if (any(fit$family$family == c("poisson", 
            "binomial"))) 
            1
        else if (fit$df.residual > 0) {
            if (any(fit$weights == 0)) 
                warning("observations with zero weight ", "not used for calculating dispersion")
            sum(fit$weights * fit$residuals^2)/fit$df.residual
        }
        else Inf
        p <- fit$rank
        if (p > 0) {
            p1 <- 1:p
            Qr <- fit$qr
            aliased <- is.na(coef(fit))
            coef.p <- fit$coefficients[Qr$pivot[p1]]
            covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            dimnames(covmat.unscaled) <- list(names(coef.p), 
                names(coef.p))
            fit$var <- dispersion * covmat.unscaled
        }
    }
# back-transformation
param <- fp.rescale(fit)
fit$coefficients <- param$coefficients
fit$var <- param$var
# Wald test (Coxph)
    if (cox) {
        if (length(fit$coefficients)) {
           if (is.null(fit$wald.test)) {
             nabeta <- !is.na(fit$coefficients)
             if (is.null(init)) 
                temp <- fit$coefficients[nabeta]
             else temp <- (fit$coefficients - init)[nabeta]
            }
            if (exists("coxph.wtest")) 
               tester <- get("coxph.wtest")
            else tester <- getFromNamespace("coxph.wtest", "survival")
            fit$wald.test <- tester(fit$var[nabeta, nabeta], temp, 
               .Machine$double.eps^0.75)$test
		}
    }
    if (x) 
        fit$X <- X
    if (y) 
        fit$y <- Y
    fit$terms <- Terms
lhs <- strsplit(as.character(fit$terms), "~")[[2]]
#varc <- attr(fit$terms,"variables")[-1]
pstvars <- NULL
if(cox) {
 stvars <- as.character(attr(fit$terms,"variables"))[-1]
 pstvars <- grep("strata\\(", stvars)
}
vars <- xnames
tvars1 <- dimnames(fit$fptable[fit$fptable$df.final!=0,])[[1]]
if(length(tvars1)) {      # are some vars selected?
 spos <- match(tvars1,dimnames(fit$scale)[[1]])
 scale <- fit$scale[spos,,drop=FALSE]
 fpos <- match(tvars1,dimnames(fit$fptable)[[1]])
 fptable <- fit$fptable[fpos,,drop=FALSE]
 tvars <- ifelse(scale[,1]!=0, paste("(",tvars1,"+",scale[,1],")",sep=""), tvars1)
 for(iv in seq(tvars)) {     # if so, then work along 
 tv <- match(tvars1[iv], vars[fp.pos])
  if(length(tv)) if(!is.na(tv))  # are there fps selected?
  {
    shift <- scale[,1][iv]
    if(fptable$power2[iv]!=".") 
	{
	  if(fptable$power1[iv]!=0) 
	  {
		if(fptable$power2[iv]==as.character(fptable$power1[iv]) | fptable$power2[iv]==0) 
		   tvars[iv] <- paste("I(",tvars[iv], "^",fptable$power1[iv], ")+log(",tvars[iv],")",sep="", collapse="")
	    else 
		   tvars[iv] <- paste("I(",tvars[iv], "^",fptable$power1[iv], ")+I(",tvars[iv], "^",fptable$power2[iv],")",sep="", collapse="")
      } else {
		if(fptable$power2[iv]==0) 
		   tvars[iv] <- paste("log(",tvars[iv],")",sep="", collapse="")
	    else 
		   tvars[iv] <- paste("log(",tvars[iv], ")+I(",tvars[iv], "^",fptable$power2[iv],")",sep="", collapse="")
	  }
	} else {
	 if(fptable$power1[iv]!=0) tvars[iv] <- paste("I(",tvars[iv], "^",fptable$power1[iv], ")",sep="", collapse="")
	 else tvars[iv] <- paste("log(",tvars[iv],")",sep="", collapse="")
    }
   } else {   # non-fp vars?
    if(length(vars[-fp.pos])) {
	 tv <- pmatch(vars[-fp.pos], tvars[iv])
     tvars[iv] <- vars[-fp.pos][which(!is.na(tv))]
	} 
  }
}
if(length(pstvars)) # cox - strata?
 rhs <- paste(c(stvars[pstvars],tvars), collapse="+")
else
 rhs <- paste(tvars, collapse="+")
formula <- eval(parse(text=paste(lhs,"~",rhs)))
} else {
formula <- eval(parse(text=paste(lhs,"~ 1")))
}
fit$formula <- formula
fit$call <- call
if(cox) {
 fit$fit <- coxph(formula=formula, data=data, subset = subset, na.action = na.action, init = init, method = method, x = x, y = y)
} else {
 fit$fit <- glm(formula = formula, family = family, data = data, subset = subset, na.action = na.action, x = x, y = y)
 fit$qr <- fit$fit$qr; fit$R <- fit$fit$R;  fit$effects <- fit$fit$effects
}
fit
}
