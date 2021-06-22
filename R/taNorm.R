#' @title Total Area Normalisation
#' @description Total area normalisation (TA) scales spectra so that they all have a total integral of one
#' @details
#' ### How It Works:
#' 1. The integral of each experimental spectra in X is calculated
#' 2. Each intensity of a spectrum is divided by it's integral which reduces it's total area to one
#' ### Background Information:
#' * TA was once considered the gold standard in NMR spectra normalisation.
#' * TA is also referred to as constant sum or integral normalisation
#' * TA is useful for spectra that have similar peaks.
#' ### Advantages:
#' * TA is a straightforward normalisation method with limited computational requirements making it quick and uncomplicated.
#' * Combining both TA and another normalisation method could also provide more robust results.
#' ### Limitations:
#' * Large variation in the peaks present in a spectra severly impact TA's performance.
#' * For instance, some spectra have large peaks created by drug metabolites, etc., which increase the total area of the spectrum.
#' * This large total area will then dramatically down scale spectra and give misleading results.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @return A list of
#' 1. The normalised X matrix in the first list element, and
#' 2. A numerical array of the corresponding dilution factors.
#' * Following the example below will extract the results quickly and easily.
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#' @examples
#' Xt <- taNorm(X)
#' Xta <- Xt[[1]]
#' taDilf <- Xt[[2]]
#' @export

taNorm <- function(X){
  Xa <- apply(X, 1, sum)
  Xta <- t(sapply(1:nrow(X), function(x){
    X[x, ]/Xa[x]
  }))
  return(list(Xa, Xta))
}
