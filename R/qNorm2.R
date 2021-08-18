#' @title Quantile Normalisation 2
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

qNorm2 <- function(X, uv_used = 'mode'){
  if (is.null(uv_used)){
    stop("uv_used was left blank. Please specify which univariate variable you want to use to calculate the dilf. 'mode' and 'median' accepted.")
  }
  if (!uv_used == 'mode' & !uv_used == 'median'){
    stop("uv_used not specified correctly. Only 'mode' and 'median' are accepted.")
  }
  cat('\033[0;34mCalculating Quantile Means... \033[0m')
  Xs <- t(apply(X, 1, sort))
  Xr <- t(apply(X, 1, rank))
  Xm <- apply(Xs, 2, mean)
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mReassigning Intensity Values... \033[0m')
  Xq <- t(sapply(1:nrow(Xr), function(i){
    Xm[Xr[i,]]
  }))
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mCalculating Quotients... \033[0m')
  q <- t(sapply(1:nrow(X), function(j){
      X[j,]/Xq[j,]
  }))
  if (uv_used == "mode"){
    cat('\033[0;34mUsing the Mode... \033[0m')
    dilf <- sapply(1:nrow(q), function(y){
      i <- q[y,]
      den <- suppressWarnings(density(log10(i)[!is.nan(log10(i)) & !is.na(log10(i)) & !is.infinite(log10(i))]))
      dilf <- 10^(den$x[which.max(den$y)])
      return(dilf)
    })
  }
  if (uv_used == "median"){
    cat('\033[0;34mUsing the Median... \033[0m')
    dilf <- sapply(1:nrow(q), function(z){
      i <- q[z]
      dilf <- median(i)
      return(dilf)
    })
  }
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mNormalising X... \033[0m')
  Xn <- t(sapply(1:nrow(X), function(a){
    X[a,]/dilf[a]
  }))
  rownames(Xn) <- rownames(X)
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_q2", Xn, envir = .GlobalEnv)
  assign("dilf_q2", dilf, envir = .GlobalEnv)
}
