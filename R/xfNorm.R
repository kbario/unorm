#' @title Normalisation based on an External Factor
#' @param X a numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param xfactor a numerical array of values that correspond to the dilution coefficients of the spectra
#' @return a list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @details https://doi.org/10.1021/ac051632c
#' @author \email{kylebario1@@gmail.com}
#' @examples
#'   Xn <- xfNorm(X, osmo)
#'   xfDilf <- Xn[[2]]
#'   Xn <- Xn[[1]]
#' @export

xfNorm <- function(X, xfactor){
  Xn <- sapply(1:nrow(X), function(x){
    X[x, ]/xfactor[x]
  })
  dilf <- xfactor
  return(list(Xn, dilf))
}
