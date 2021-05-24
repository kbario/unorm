#' @title Total Area Normalisation
#' @param X spectral matrix with rows being the spectra and columns being the chemical shift variables
#' @return list of normalised X matrix and an array of the corresponding dilution factors
#' @details also called Constant Sum normalisation
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
