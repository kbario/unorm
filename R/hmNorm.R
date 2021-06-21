#' @title Histogram Matching
#' @description [hmNorm()] is a reference based normalisation that aims to match the histogram created from a experimental spectra with that made from a reference spectra.
#' @param X A numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @param ppm A numerical array holding the chemical shift values of the X matrix
#' @return A list with the normalised X matrix and an array of the corresponding dilution factors.
#' @details #Limitations:
#' * Histogram matching (HM) overcomes issues like peak shift that impact methods such as PQN, however, HM has limitations caused by differently shaped histograms.
#' * Because not all histograms will have the same level of skewness or kurtosis, not all matchings will be ideal which limits this methods applicability.
#' @seealso The methods paper for HM can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @examples
#' \dontrun{
#'  histMod <- hmNorm(X, ppm)
#'  hmDilf <- histMod[[2]]
#'  Xn <- histMod[[1]]}
#' @export
#' @family {Reference-Based}
#' @author \email{kylebario1@@gmail.com}
#' @importFrom stats sd
#' @importFrom metabom8 get_idx
#' @importFrom graphics hist

hmNorm <- function(X, ppm, noi = c(10,11), intensity_binwidth = 0.1, alpha_from = 0.5, alpha_to = 1.5, alpha_n = 101, use_median = F){
  std <- apply(X[,get_idx(noi, ppm)], 1, function(i){
    sd(i)
  })
  alpha <- seq(from = alpha_from, to = alpha_to, length.out = alpha_n)
  Xstar <- t(sapply(1:nrow(X), function(j){
    out = X[j, ]
    out[(X[j,] <= std[j]*5)] = NA
    return(out)
  }))
  Xstar_small <- Xstar*alpha[1]
  Zsmall <- log2(abs(Xstar_small)+1)
  Xstar_big <- Xstar*alpha[alpha_n]
  Zbig <- log2(abs(Xstar_big)+1)
  width <- ((max(Zbig, na.rm = T)-min(Zsmall, na.rm = T))/intensity_binwidth)
  br <- seq(from = min(Zsmall, na.rm = T), to = max(Zbig, na.rm = T), length.out = width)
  if (use_median){
    Xstar_median <- apply(Xstar, 2, median)
    Zmedian <- log2(abs(Xstar_median)+1)
    Hmedian <- hist(Zmedian, breaks = br, plot = F)
    matching <- sapply(1:nrow(X), function(i){
      sapply(1:length(alpha), function(j){
        star_new <- Xstar[i,]*alpha[j]
        z <- log2(abs(star_new)+1)
        H <- hist(z, breaks = br, plot = F)
        ss <- sum(((Hmedian$counts-H$counts)^2))
        return(ss)
      })
    })
    matching <- t(matching)
  } else {
    Zt <- log2(abs(Xstar[1,])+1)
    Htarget <- hist(Zt, breaks = br, plot = F)
    matching <- sapply(1:nrow(X), function(i){
      sapply(1:length(alpha), function(j){
        star_new <- Xstar[i,]*alpha[j]
        z <- log2(abs(star_new)+1)
        H <- hist(z, breaks = br, plot = F)
        ss <- sum(((Htarget$counts-H$counts)^2))
        return(ss)
      })
    })
    matching <- t(matching)
  }
  alpha_min <- sapply(1:nrow(X), function(k){
    (alpha[which.min(matching[k,])])
  })
  Xn <- t(sapply(1:nrow(X), function(l){
    X[l,]/alpha_min[l]
  }))
  return(list(Xn,alpha_min))
}
