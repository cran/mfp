fp.scale <- function(x, scaling = TRUE)
{
#
# Version 1.3     27.03.2005
#
scale <- 1; shift <- 0
  if(scaling) {
      if(min(x) <= 0) {
#        z <- sort(x)[-1] - sort(x)[ - length(x)]
        z <- diff(sort(x))
        shift <- min(z[z > 0]) - min(x)
        shift <- ceiling(shift*10)/10
      }
      else shift <- 0
    range <- max(x) - min(x)
    range <- median(x+shift)
#    scale <- 10^(sign(log10(range)) * trunc(abs(log10(range))))
    scale <- 10^(sign(log10(range)) * round(abs(log10(range))))
  }
    return(list(shift=shift, scale=scale))
}
