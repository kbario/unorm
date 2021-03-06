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
#' @param noi Takes an array that is row matched to the X matrix you are normalising with the values equaling the maximum noise estimation for each spectra respectively.
#' @return A list of
#' 1. The normalised X matrix in the first list element, and
#' 2. A numerical array of the corresponding dilution factors.
#' * Following the example below will extract the results quickly and easily.
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#' @examples
#' data(X, noi)
#' taNorm(X, noi)
#' cat(dilf_ta)
#' @export

taNorm <- function(X, noi){
  cat('\033[0;34mCalculating Dilfs... \033[0m')
  Xs <- t(sapply(1:nrow(X), function(i){
    x <- X[i,]
    x[x<noi[i]]=0
    return(x)
  }))
  Xa <- unname(apply(Xs, 1, sum))
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mNormalising X... \033[0m')
  Xta <- t(sapply(1:nrow(Xs), function(x){
    X[x, ]/Xa[x]
  }))
  rownames(Xta) <- rownames(X)
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_ta", Xta, envir = .GlobalEnv)
  assign("dilf_ta", Xa, envir = .GlobalEnv)
}
