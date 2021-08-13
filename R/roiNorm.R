#' @title Region of Interest Normalisation
#' @description This function normalises the spectra based on a specific area of the spectra. This is helpful when one region is constant across all spectra.
#' @details
#' ### How It Works:
#' 1. This function takes the values of the shift argument and finds the integral (or area) of the peak in that chemical shift region in an experimental spectra.
#' 2. Each intensity in that experimental spectra is then divided by that integral to scale it.
#' This method works much like how [taNorm()] and [creNorm()] do, except instead of using an entire spectrum's integral or the creatinine peak specifically, you are able to use what ever ppm region you wish.
#' ### Advantages:
#' * Region of Interest Normalisation (ROI) gives the data analyst more control over what region will be used which can be advantageous.
#' ### Limitations:
#' * ROI is also subject to the flaws that impact creatinine normalisation ([creNorm()]) in that, whatever factors that create variation in the region one chooses, will also impact the entire spectrum's state.
#' * Reference-based normalisations may provide more robust normalisation.
#' @family {Attribute-Based}
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param ppm A numerical array holding the chemical shift values of the X matrix
#' @param shift The numerical values defining the lower and upper regions of the Region of Interest. default = c(3,3.1).
#' @return This function returns a list with:
#' 1. The normalised X matrix, and
#' 2. A numerical array of the corresponding dilution factors.
#' * Following the example below will extract the results quickly and easily.
#' @seealso A description of Region of Interest Normalisation can be found in this paper: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' roy <- roiNorm(X, ppm, shift = c(2,2.3))
#' Xn <- roy$Xn
#' dilf <- roy$dilf
#' @importFrom metabom8 get_idx
#' @export

roiNorm <- function(X, shift = c(2.5,2.75)){
  p <- as.numeric(colnames(X))
  i <- get_idx(shift, p)
  cat('\033[0;34mCalculating Dilfs... \033[0m')
  dilf <- sapply(1:nrow(X), function(y){
    (sum(X[y,i]))
  })
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mNormalising X... \033[0m')
  Xn <- t(sapply(1:nrow(X), function(x){
    (X[x,])/(sum(X[x,i]))
  }))
  rownames(Xn) <- rownames(X)
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_roi", Xn, envir = .GlobalEnv)
  assign("dilf_roi", dilf, envir = .GlobalEnv)
}
