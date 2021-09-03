#' @title Quantile Normalisation
#' @description A method of NMR spectral normalisation where the maximum intensities
#' @details The intensities of each spectrum are ordered smallest to largest, the means of these 'quantiles' are calculated and these means are then reassigned to the ppm they match to (i.e., the ppm which contained the highest intensity will be assigned the mean of the highest values and so on)
#' @family {Attribute-Based}
#' @param X A numerical matrix containing the NMR spectra to be normalised. Rows should be the spectra and columns being the chemical shift variables
#' @return This function returns a list with:
#' 1. The normalised X matrix based on the calculated dilf in the first element,
#' 2. The quantile normalised X matrix in the second element, and
#' 3. The dilf calculated in the third element.
#' @seealso The methods paper that describes quantile normalisation: \url{https://doi.org/10.1007/s11306-011-0350-z}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' data(X)
#' qNorm(X)
#' @export

qNorm <- function(X){
  cat('\033[0;34mCalculating Quantile Means... \033[0m')
  Xs <- t(apply(X, 1, sort))
  Xr <- t(apply(X, 1, rank))
  Xm <- apply(Xs, 2, mean)
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mReassigning Intensity Values... \033[0m')
  Xq <- t(sapply(1:nrow(Xr), function(i){
    Xm[Xr[i,]]
  }))
  rownames(Xq) <- rownames(X)
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_q", Xq, envir = .GlobalEnv)
}
