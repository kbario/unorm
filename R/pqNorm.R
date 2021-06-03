#' @title Probabilistic Quotient Normalisation
#' @param X spectral matrix with rows being the spectra and columns being the chemical shift variables
#' @return list of normalised X matrix and an array of the corresponding dilution factors
#' @details https://doi.org/10.1021/ac051632c
#' @author kylebario1@@gmail.com
#' @example
#' Xn <- pqNorm(X)
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
