\name{plot.mfp}
\alias{plot.mfp}

\title{Plots for mfp objects}

\description{This function draws two plots: (i) the linear predictor function and (ii) partial residuals together with a lowess smooth. 
For Cox models also smoothed martingale based residuals of the null 
model are plotted against the predictor.
}
   
\usage{plot.mfp(x, ...)}

\arguments{
 \item{x}{
	object representing a fitted \code{mfp} model.}
 \item{...}{ 
        further arguments.}
}

\value{}

\examples{
\dontrun{
	data(GBSG)
	f <- mfp(Surv(rfst, cens) ~ fp(age, df = 4, select = 0.05)
                 + fp(prm, df = 4, select = 0.05), family = cox, data = GBSG)
	plot(f)
}
}

\keyword{ }
