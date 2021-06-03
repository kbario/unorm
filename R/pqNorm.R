#' @title Probabilistic Quotient Normalisation
#' @param X A matrix containing the NMR spectra with rows contain all data points from one spectra and the columns containing the chemical shift variables.
#' @return A list with the normalised X matrix and an array of the corresponding dilution factors.
#' @details The methods paper can be found here: https://doi.org/10.1021/ac051632c
#' @author \email{kylebario1@@gmail.com}
#' @example
#'     Xn <- pqNorm(X)
#' @export
#' @importFrom stats median density

pqNorm <- function(X){
  Xm <- apply(X, 2, median)
  dilf <- sapply(1:nrow(X), function(x){
    d <- density((X[x,]/Xm), from = -10, to = 10)
    dilf <- d$x[which.max(d$y)]
    return(dilf)
  })
  Xn <- sapply(1:length(dilf), function(x){
    X[x,]/dilf[x]
  })
  return(list(Xn, dilf))
}
