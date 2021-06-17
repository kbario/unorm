#' @title Region of Interest Normalisation (ROI)
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param ppm A numerical array holding the chemical shift values of the X matrix
#' @param shift The numerical values defining the lower and upper regions of the Region of Interest. default = c(3,3.1).
#' @return A list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @example
#'   Xn <- roiNorm(X, ppm, shift = c(2,2.3))
#'   roiDilf <- Xn[[2]]
#'   Xn <- Xn[[1]]
#' @details http://dx.doi.org/10.1007/s11306-018-1400-6
#' @author \email{kylebario1@@gmail.com}
#' @export
#' @importFrom metabom8 get_idx

roiNorm <- function(X, ppm, shift = c(3,3.1)){
  i <- get_idx(range = shift, ppm)
  a <- apply(X[, i], 1, sum)
  Xn <- sapply(1:nrow(X), function(n){
    X[n,]/a
    return(Xn)
  })
  return(list(Xn, a))
}
