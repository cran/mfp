## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## -----------------------------------------------------------------------------
library(mfp)

str(mfp)

## -----------------------------------------------------------------------------
str(fp)

## -----------------------------------------------------------------------------
data(GBSG)
str(GBSG)

## -----------------------------------------------------------------------------
f <- mfp(Surv(rfst, cens) ~ strata(htreat)+age+fp(tumsize)+fp(posnodal)+fp(prm)+fp(esm)
        +menostat+tumgrad, family = cox, data = GBSG, select=0.05, verbose=TRUE)

## -----------------------------------------------------------------------------
summary(f)

## -----------------------------------------------------------------------------
f$fptable

## ----fig2---------------------------------------------------------------------
vizmfp <- predict(f, type = "terms", terms = "posnodal", seq = list(1:50), ref = list(5))

plot(vizmfp$posnodal$variable, exp(vizmfp$posnodal$contrast), type = "n", log = "y",
     xlab = "posnodal", ylab = "Hazard Ratio", ylim = c(0.1, 5))
polygon(x = c(vizmfp$posnodal$variable, rev(vizmfp$posnodal$variable)),
        y = exp(c(vizmfp$posnodal$contrast - 1.96 * vizmfp$posnodal$stderr,
              rev(vizmfp$posnodal$contrast + 1.96 * vizmfp$posnodal$stderr))),
        col = "grey", border = NA)
grid()
lines(vizmfp$posnodal$variable, exp(vizmfp$posnodal$contrast), type = "l", col = 4, lwd = 2)


## ----fig1---------------------------------------------------------------------
pf <- survfit(f$fit)  
plot(pf, col=c("red","green"), xlab="Time (years)", ylab="Recurrence free survival rate", xscale=365.25)

