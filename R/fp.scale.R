fp.scale <- function(x, scaling)
{
#
#   Version 1.0.0   19 Oct 1999
#
scale <- 1; shift <- 0
  if(scaling) {
      if(min(x) <= 0) {
        z <- sort(x)[-1] - sort(x)[ - length(x)]
        shift <- min(z[z > 0]) - min(x)
    }
    else shift <- 0
    range <- max(x) - min(x)
    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
  }
    return(list(shift=shift, scale=scale))
}
