### R code from vignette source 'mfp.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt=">", width=80)
set.seed(20040804)


###################################################
### code chunk number 2: mfp.Rnw:94-95
###################################################
library(mfp)


###################################################
### code chunk number 3: mfp.Rnw:105-106
###################################################
str(mfp)


###################################################
### code chunk number 4: mfp.Rnw:128-129
###################################################
str(fp)


###################################################
### code chunk number 5: mfp.Rnw:167-169
###################################################
data(GBSG)
str(GBSG)


###################################################
### code chunk number 6: mfp.Rnw:194-196
###################################################
f <- mfp(Surv(rfst, cens) ~ strata(htreat)+age+fp(tumsize)+fp(posnodal)+fp(prm)+fp(esm)
        +menostat+tumgrad, family = cox, data = GBSG, select=0.05, verbose=TRUE)


###################################################
### code chunk number 7: mfp.Rnw:206-207
###################################################
summary(f)


###################################################
### code chunk number 8: mfp.Rnw:212-213
###################################################
f$fptable


###################################################
### code chunk number 9: FIG1
###################################################
pf <- survfit(f$fit)  
plot(pf, col=c("red","green"), xlab="Time (years)", ylab="Recurrence free survival rate", xscale=365.25)


