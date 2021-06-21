#' @title Total Area Normalisation
#' @description [taNorm()] finds the total integral of the spectra and divides each spectra by its integral.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @return A list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @details
#' Total area normalisation (TA) was once considered the gold standard in NMR spectra normalisation.
#' TA is useful for spectra that have similar peaks. If spectra differ a lot from each other - for instance, some have large peaks created by drug metabolites, etc. - then the total area will not be a suitable normalisation method.
#' TA is also referred to as constant sum or integral normalisation
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#' @examples
#' \dontrun{
#'     Xn <- taNorm(X)
#'     taDilf <- Xn[[2]]
#'     Xn <- Xn[[1]]}
#' @export

taNorm <- function(X){
  Xa <- apply(X, 1, sum)
  Xta <- t(sapply(1:nrow(X), function(x){
    X[x, ]/Xa[x]
  }))
  return(list(Xa, Xta))
}
