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
#' @param ppm A numerical array holding the chemical shift values of the X matrix. Should be column matched to X provide.
#' @param noi The concatenated ppm values defining the lower and upper chemical shift values of the noise. There should be **no** signal between these chemical shift values. Also be sure to have not cut out these regions of noise in the [preprocessing()] function.
#' @param intensity_binwidth This argument dictates the width of the bins. The average span of intensities is from `10`-`30`, meaning that an `intensity_binwidth` of `0.1` would give you 200 bins.
#' @param alpha_from `alpha` is the scaling factor that the experiment histograms will be scaled with. `alpha_from` defines what the lowest value of `alpha` will be. Can be as low as `0.1` but `0.5` is default. Consider how much dilution variation there is in your samples before setting `alpha_from` extremely low, the more `alpha`s to test, the longer the computation time will be.
#' @param alpha_to This defines the largest value of `alpha`. Default is `1.5.` As above, consider the dilution variation of the samples. The larger `alpha_to` is the longer the computation time.
#' @param alpha_n The number of `alpha`s you wish to test. `alpha_n`s default is `101` which will give you `alpha`s incrementing by `0.5.` The larger `alpha_n`, the longer the computing will take and more computing power you will need.
#' @param use_median This argument dictates whether the function will calculate the median and use that as the reference spectrum or not. If set to `FALSE`, the first sample will be used as the reference. This practice is outlined in the methods paper (see 'See also')
#' @return A list with:
#' 1. The normalised X matrix in the first list element, and
#' 2. An array of the corresponding dilution factors.
#' * Following the example below will extract the results quickly and easily.
#' @seealso
#' * The methods paper for HM can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' * The paper discussing HM limitations with noise can be found here: \url{http://dx.doi.org/10.1007/s11306-018-1400-6}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#'  Xh <- hmNorm(X, ppm)
#'  Xhm <- Xh[[1]]
#'  hmDilf <- Xh[[2]]
#' @importFrom stats sd
#' @importFrom metabom8 get_idx
#' @importFrom graphics hist
#' @export

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
