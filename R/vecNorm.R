#' Vector Length Normalisation
#' @description This function aims to normalise NMR spectra by their vector length.
#' @details This function operates by squaring all values of the spectrum, summing them, and square rooting the result. Then each element from that spectrum is divided by this vector length.
#'
#' @param X The spectr(a/um) to be normalised
#'
#' @return The return value of `vecNorm()` will be of same shape as the X argument parsed to it, either an array or matrix of normalised X
#' @export
#'
#' @author \email{kylebario1@@gmail.com}
#' @family {Attribute-Based}
#'
#' @seealso The methods paper that describes this method of normalisation can be found here: \url{https://doi.org/10.1021/ac051632c}
#'
#' @examples
#' data(X)
#' vecNorm(X)
#' cat(dilf_vec)
vecNorm <- function(X){
  if (is.null(nrow(X))){
    cat('\033[0;34mCalculating Dilf... \033[0m')
    l <- sqrt(sum(X^2))
    cat('\033[1;32mDone.\n\033[0m')
    cat('\033[0;34mNormalising X... \033[0m')
    Xn <- X/l
    cat('\033[1;32mDone.\n\033[0m')
  } else {
    cat('\033[0;34mCalculating Dilfs... \033[0m')
    l <- apply(X, 1, function(x){
      sqrt(sum(x^2))
    })
    cat('\033[1;32mDone.\n\033[0m')
    cat('\033[0;34mNormalising X... \033[0m')
    Xn <- sapply(1:length(l), function(y){
      X[y,]/l[y]
    })
    cat('\033[1;32mDone.\n\033[0m')
  }
  assign("X_vec", Xn, envir = .GlobalEnv)
  assign("dilf_vec", l, envir = .GlobalEnv)
}
