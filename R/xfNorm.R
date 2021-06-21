#' @title Normalisation based on an External Factor
#' @description One of the best ways to normalise urine spectra is by obtaining information about the dilution of the sample from an external source such as freezing point osmilality. The result from the external source can be used to scale spectra and give an over-all accurate result.
#' @details Normalisation based on an external factor is a more robust method of normalisation for urine spectra especially using the gold standard of dilution measurement - freezing-point depression osmometry - but this requires more effort and time compared to other normalisation methods like PQN.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param xfactor A numerical array of values that correspond to the dilution coefficients of the spectra
#' @return This function returns a list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Reference-Based}
#' @examples
#' \dontrun{
#'   Xn <- xfNorm(X, osmo)
#'   xfDilf <- Xn[[2]]
#'   Xn <- Xn[[1]]}
#' @export

xfNorm <- function(X, xfactor){
  Xn <- sapply(1:nrow(X), function(x){
    X[x, ]/xfactor[x]
  })
  dilf <- xfactor
  return(list(Xn, dilf))
}
