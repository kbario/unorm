#' Baseline Correction
#'
#' @param X The X matrix you wish to correct the baseline of
#' @param lambda The smoothing parameter
#' @param it_max The maximum iterations before bl stops
#'
#' @return A matrix containing spectra with corrected baselines
#' @export
#'
#' @family {Data_Manipulation}
#'
#' @examples
#' path = system.file('extdata', package = 'unorm')
#' X_bl <- bl(X)
bl <- function(X, lambda = 1e+07, it_max = 30)
{
  if (any(is.na(X))) {
    cat("\033[0;33mReplacing NA's in X with zeros... \033[0m")
    X[is.na(X)]=0
  }
  if (is.null(ncol(X))){
    X <- t(X)
  }
  Xb <- t(apply(X, 1, function(x) {
    x - ptw::asysm(x, maxit = it_max, lambda = lambda)
  }))
  return(Xb)
}
