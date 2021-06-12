#' @title Quantitle Normalisation
#' @param X A matrix containing the NMR spectra with rows contain all data points from one spectra and the columns containing the chemical shift variables.
#' @return A list with the normalised X matrix and an array of the corresponding dilution factors.
#' @details
#' @author \email{kylebario1@@gmail.com}
#' @export

qNorm <- function(X){
  Xs <- t(apply(X, 1, sort))
  Xo <- t(apply(X, 1, order, decreasing = T))
  Xmean <- apply(Xs, 2, mean)
  Xn <- t(apply(Xo, 1, function(i, imean = Xmean){#browser()
    imean[i]
  }))
  return(Xn)
}

# Xb <- binning(X, ppm, width = 0.05)
# X <- Xb
#
# ppb=as.numeric(colnames(Xb))
# matspec(Xn, ppb, interactive = F)
#
#
#
# plot(Xb[rownames(Xb)=='10', ]~ppb, type='l')
# points(Xn[rownames(Xb)=='10', ]~ppb, type='l', col='blue')
#
#
#
#
# rm(list = ls())
