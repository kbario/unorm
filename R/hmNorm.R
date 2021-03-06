#' @title Histogram Matching
#' @description Histogram Matching (HM) is a reference based normalisation that aims to match the histogram created from a experimental spectra with that made from a reference spectra.
#' @details
#' ### How It Works:
#' 1. A more complex normalisation method compared to others in this package, HM creates a histogram of the intensities of an experimental spectrum.
#' 2. This area of this histogram subtracted from that of a reference histogram, calculated from a median spectrum.
#' 3. The experimental spectra is scaled by a factor `alpha` which alters the intensities and a new histogram is made and step 2. is repeated.
#' 4. The aim is to minimise the difference between the reference histogram and the experiment histogram and the alpha value that acheives this minimum is the dilution coefficient.
#' 5. This is iterated over for all experimetal spectra
#' ### Advantages:
#' * HM is not affected by peak shift like PQN ([pqNorm()]) is because it does not take into account the chemical shift of the intensity, on the value of the intensity.
#' * HM also is not as impaacted by noise because noise is removed before histograms are produced (this is not touched on in the above section)
#' ### Limitations:
#' * Not all histograms are created equal, some will have more kurtosis or skewness which limits the degree to which differences can be minimised.
#' * Studies have also found that HM does not handle noise well and does not recover signal like many other normalisation methods do. (see 'See also')
#' @family {Reference-Based}
#' @param X A numerical matrix containing the NMR spectra to be normalised. Rows should be the spectra and columns being the chemical shift variables
#' @param noi The array of maximum noise estimations produced from the function [noise()] or provided by the user. Must match X rows.
#' @param int_binwid This argument dictates the width of the bins. The average span of intensities is from `10`-`30`, meaning that an `intensity_binwidth` of `0.1` would give you 200 bins.
#' @param use_median This argument dictates whether the function will calculate the median and use that as the reference spectrum or not. If set to `FALSE`, the first sample will be used as the reference. This practice is outlined in the methods paper (see 'See also')
#' @param alpha The lower and upper bounds that the golden selection search will search between
#' @param tol This defines the tolerance or level of precision the golden selection search will search until. (i.e., it will search until the bounds are `tol` distance apart)
#' @return A list with:
#' 1. The normalised X matrix in the first list element, and
#' 2. An array of the corresponding dilution factors.
#' * Following the example below will extract the results quickly and easily.
#' @seealso
#' * The methods paper for HM can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' * The paper discussing HM limitations with noise can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' # it is mandatory to fill args X, noi, and use_median for `hmNorm()` to work.
#' data(X, noi)
#' hmNorm(X, noi, use_median = TRUE, alpha = c(0.1, 2.5))
#' cat(dilf_hm)
#' @importFrom stats sd
#' @importFrom graphics hist
#' @importFrom stats median
#' @export

hmNorm <- function(X, noi, int_binwid = 0.1, alpha = c(0.1, 2), use_median = NULL, tol = 1e-5){
  if (is.null(use_median)){
    stop("Please define which reference to use. To use median: use_median = T; to use first spectra: use_median = F.")
  }
  cat('\033[0;34mPrepping the spectra\n\033[0m')
  Xs <- t(sapply(1:nrow(X), function(j){
    Xc <- X[j,]
    n <- noi[j]
    Xc[(Xc <= n)] = NA
    return(Xc)
  }))
  Zs <- log2((Xs*alpha[1])+1)
  Zb <- log2((Xs*alpha[2])+1)
  br <- seq(from = min(Zs, na.rm = T), to = max(Zb, na.rm = T), length.out = ((max(Zb, na.rm = T)-min(Zs, na.rm = T))/int_binwid))
  if (use_median){
    cat('\033[0;34mUsing the median spectra as reference\n\033[0m')
    cat('\033[0;34mCalculating the dilfs... \033[0m')
    Xsm <- apply(Xs, 2, median)
    Zm <- log2(Xsm+1)
    m <- table(cut(Zm, breaks = br))
    dilfs <- sapply(1:nrow(Xs), function(i){
      a <- alpha[1]
      b <- alpha[2]
      gr <- (sqrt(5)+1)/2
      c <- b-(b-a)/gr
      d <- a+(b-a)/gr
      while(abs(b-a)>tol){
        cs <- table(cut(log2(Xs[i,]*c), breaks = br))
        fc <- sum((m-cs)^2)
        ds <- table(cut(log2(Xs[i,]*d), breaks = br))
        fd <- sum((m-ds)^2)
        if (fc<fd){
          b <- d
        } else {
          a <- c
        }
        c <- b-(b-a)/gr
        d <- a+(b-a)/gr
      }
      return((b+a)/2)
    })
    cat('\033[1;32mDone.\n\033[0m')
  } else {
    cat('\033[0;34mUsing the first spectra as reference\n\033[0m')
    cat('\033[0;34mCalculating dilfs... \033[0m')
    m <- table(cut((log2(Xs[1,]+1)), breaks = br))
    dilfs <- sapply(1:(nrow(Xs)), function(i){
      if (i==1){
        return(i)
      } else {
        a <- alpha[1]
        b <- alpha[2]
        gr <- (sqrt(5)+1)/2
        c <- b-(b-a)/gr
        d <- a+(b-a)/gr
        while(abs(b-a)>tol){
          cs <- table(cut(log2(Xs[i,]*c), breaks = br))
          fc <- sum((m-cs)^2)
          ds <- table(cut(log2(Xs[i,]*d), breaks = br))
          fd <- sum((m-ds)^2)
          if (fc<fd){
            b <- d
          } else {
            a <- c
          }
          c <- b-(b-a)/gr
          d <- a+(b-a)/gr
        }
      }
      return((b+a)/2)
    })
    cat('\033[1;32mDone.\n\033[0m')
  }
  cat('\033[0;34mNormalising X... \033[0m')
  Xn <- t(sapply(1:nrow(X), function(l){
    X[l,]/dilfs[l]
  }))
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_hm", Xn, envir = .GlobalEnv)
  assign("dilf_hm", dilfs, envir = .GlobalEnv)
}
