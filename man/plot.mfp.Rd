\name{plot.mfp}
\alias{plot.mfp}

\title{Plots for mfp objects}

\description{This function draws two plots: (i) the linear predictor function and (ii) partial residuals together with a lowess smooth. 
For Cox models also smoothed martingale based residuals of the null 
model are plotted against the predictor.
}
   
\usage{plot.mfp(x, var=NULL, ref.zero=TRUE, ...)}

\arguments{
 \item{x}{
	object representing a fitted \code{mfp} model.}
 \item{var}{
	the variable for which plots are desired. 
      By default, plots are produced in turn for each variable of a model.}
 \item{ref.zero}{
        subtract a constant from X beta before plotting so that the
          reference value of the 'x'-variable yields 'y=0'.}
 \item{...}{ 
        further arguments.}
}

\value{}

\examples{
	data(GBSG)
	f <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05)
                 + fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)
	par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(1.5,0.75,0))
      plot(f, var="age")
}

\keyword{ }
