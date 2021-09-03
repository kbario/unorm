#' @title Normalisation based on an External Factor
#' @description This function utilises an array of dilution coefficients acquired externally to normalise the spectra.
#' @details
#' ### How It Works:
#' 1. All the intensities of a spectrum are divided by its external dilution coefficient.
#' 2. This is done to all spectra with their corresponding dilution factors
#' ### Advantages:
#' * Reliable: dilution estimations from external sources such as freezing-point depression osmometry are highly accurate, reducing the error in normalisation.
#' * Computationally, this is not as intense as other normalisation methods like [pqNorm()] and [hmNorm()] and will take limited power to perform
#' ### Limitations:
#' * This method requires having a machine to provide the external measure of dilution such as a freezing-point depression osmometer.
#' * Acquiring data on dilution with a manual freezing-point depression osmometer is a long and tedious process with each sample taking on average 5 minutes. Automatic osmometers are available but are expensive so this method is limited by logistics of actually acquiring the data.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param xfactor A numerical array of values that correspond to the dilution coefficients of the spectra. Must be row-matched to the provided X matrix.
#' @return This function returns a list with:
#' 1. The normalised X matrix in the first list element, and
#' 2. A numerical array of the corresponding dilution factors in the second.
#' * Following the example below will extract the results quickly and easily.
#' @seealso \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Reference-Based}
#' @examples
#' data(X)
#' xfNorm(X, c(300, 300))
#' @export

xfNorm <- function(X, xfactor){
  cat('\033[0;34mNormalising X... \033[0m')
  Xn <- t(sapply(1:nrow(X), function(x){
    X[x, ]/xfactor[x]
  }))
  rownames(Xn) <- rownames(X)
  cat('\033[1;32mDone.\n\033[0m')
  dilf <- xfactor
  assign("X_xf", Xn, envir = .GlobalEnv)
  assign("dilf_xf", dilf, envir = .GlobalEnv)
}
