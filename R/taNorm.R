#' @title Total Area Normalisation (TAN/CSN)
#' @param X a numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @return a list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @details also called Constant Sum normalisation https://doi.org/10.1021/ac051632c
#' @author \email{kylebario1@@gmail.com}
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
