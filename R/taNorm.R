#' @title Total Area Normalisation
#' @param X A matrix containing the NMR spectra with rows contain all data points from one spectra and the columns containing the chemical shift variables.
#' @return A list with the normalised X matrix and an array of the corresponding dilution factors.
#' @details also called Constant Sum or Integral normalisation
#' @author kylebario1@@gmail.com
#' @example
#' Xn <- taNorm(X)
#' @export

taNorm <- function(X){
  Xa <- apply(X, 1, sum)
  Xta <- t(sapply(1:nrow(X), function(x){
    X[x, ]/Xa[x]
  }))
  return(list(Xa, Xta))
}
