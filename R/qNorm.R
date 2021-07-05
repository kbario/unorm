#' @title Quantile Normalisation
#' @description A method of NMR spectral normalisation where the maximum intensities
#' @details ### How It Works:
#' ### Advantages:
#' ### Limitations:
#' @family {Attribute-Based}
#' @param X A numerical matrix containing the NMR spectra to be normalised. Rows should be the spectra and columns being the chemical shift variables
#' @return This function returns a list with:
#' 1. The normalised X matrix based on the calculated dilf in the first element,
#' 2. The quantile normalised X matrix in the second element, and
#' 3. The dilf calculated in the third element.
#' @seealso The methods paper that describes quantile normalisation: \url{https://doi.org/10.1007/s11306-011-0350-z}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' Xqn <- qNorm(X)
#' Xn <- Xqn$Xn
#' Xq <- Xqn$Xq
#' qDilf <- Xqn$dilf
#' @export

qNorm <- function(X){
  cat('\033[0;34mCalculating Xq... ')
  Xs <- t(apply(X, 1, sort))
  Xr <- t(apply(X, 1, rank))
  Xm <- apply(Xs, 2, mean)
  Xq <- t(sapply(1:nrow(X), function(x){
    Xm[Xr[x,]]
  }))
  cat('\033[1;32mDone.\n')
  Xsum <- t(apply(X, 1, sum))
  Xqsum <- t(apply(Xq, 1, sum))
  cat('\033[0;34mCalculating Dilfs... ')
  dilf <- sapply(1:nrow(X), function(h){
    Xsum[h]/Xqsum[h]
    })
  cat('\033[1;32mDone.\n')
  cat('\033[0;34mCalculating Xn... ')
  Xn <- t(sapply(1:nrow(X), function(x){
    X[x,]/dilf[x]
  }))
  cat('\033[1;32mDone.\n')
  return(list(Xn = Xn, Xq = Xq, dilf = dilf))
}

