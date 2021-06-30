#' @title Probabilistic Quotient Normalisation
#' @description PQN is currently the gold standard method used to normalise NMR spectra.
#' @details
#' ### How It Works:
#' PQN works by normalising experimental spectra in the provided X matrix in relation to a reference spectrum.
#' 1. [pqNorm()] creates the reference automatically by calculating the median spectrum of X as outlined in the initial methods paper (see 'See also')
#' 2. [pqNorm()] derives a quotient for each ppm value within the limits of shift in a experimental spectrum by dividing it's intensities with that of the reference spectrum's.
#' 3. The most frequently occurring (your choice of median or mode) quotient is calculated and is said to be the dilution coefficient (`dilf`) of that spectrum.
#' 4. The sample is then scaled with this `dilf` and will be comparable with all other spectra normalised with the reference spectrum.
#' ### Background Information:
#' * PQN is currently the gold standard for normalising NMR spectra.
#' * Multiple studies have tested its validity.
#' * PQN reliably normalises spectra and handles large amounts of noise.
#' * \strong{Using [pqNorm()] on urine spectra will give best results out of the all methods except [xfNorm()]}
#' ### Advantages:
#' * PQN is not impacted by outlining signals in the way that total area normalisation ([taNorm()]) is because it aims to find the most frequently occuring (median) `dilf`.
#' * PQN is not reliant on one signal such as creatinine normalisation ([creNorm()]) meaning that any factors that impact that one signal have less weight and don't affect the calculation of the `dilf` as much.
#' ### Limitations:
#' * PQN's biggest limitation is that it operates on the assumption that the majority of the spectra will stay constant so hugely varying spectra may be suboptimally normalised.
#' * PQN is also vulnerable to peak shift because it compares intensities of the same ppm variable. If a peak is left- or right-shifted, the apex of one signal will not be compared against the apex of another and thus produce convoluted results.
#' @family {Reference-Based}
#' @param X The numerical matrix containing the NMR data you wish to normalise. **This should be a preprocessed matrix** with baseline correction, tsp calibration and non-quantitative region removal performed on it. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @param ppm An array of chemical shift variables. ppm should be column matched to the X matrix you are normalising.
#' @param shift The concatenated ppm values that define the lower and upper bounds of the region PQN should be performed on. Regions containing noise are not useful and ideally should be omitted.
#' @param use_ta Requires a boolean `TRUE` or `FALSE` if total area normalisation should be performed on the spectra before PQN is.
#' @param uv_used PQN utilises finding the median or the mode, which are both *U*ni*v*ariate methods. Recognises either the string 'median' or 'mode' to instruct which method to use.
#' @param width Represents the bandwidth used in the 'mode' method. Only required when using `uv_used = 'mode'`. Default = 0.1. Look at the help section of [stats::density()] for more information.
#' @param noise Requires the concatenated ppm values that define the lower and upper noise regions
#' @param binning Determines if the data should be binned before performing the normalisation and takes either a `TRUE` or `FALSE` statement.
#' @return The output of this function is a list containing:
#' 1. The normalised version of X in the first element and
#' 2. A numerical array of the corresponding dilution factors calculated by the function.
#' * Following the example below will extract the results quickly and easily.
#' @seealso The methods paper first describing PQN can be found here: \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' pq <- pqNorm(X)
#' Xn <- pq$Xn
#' dilf <- pq$dilf
#' @importFrom metabom8 get_idx
#' @importFrom stats sd density median
#' @export

pqNorm <- function(X, ppm, binning = T, use_ta = F, uv_used = 'mode', width = 0.7 ,shift = c(0.5,9.5), noise = c(9.5,11)){
  if (use_ta){
    X <- t(sapply(1:nrow(X), function(x){
      (X[x,])/(sum(X[x,]))
    }))
  }
  idx <- get_idx(shift, ppm)
  Xm <- apply(X, 2, median)
  mrm <- 5*sd(Xm[get_idx(noise, ppm)])+(mean((Xm[get_idx(noise, ppm)]), trim = 0.05))
  Xm <- Xm[idx]
  Xm[Xm<=mrm] = NA
  dilf <- sapply(1:nrow(X), function(y){
    Xc <- X[y,]
    int_rm <- 5*sd(Xc[get_idx(noise, ppm)])+mean(Xc[get_idx(noise, ppm)], trim = .05)
    Xc[Xc<=int_rm]=NA
    Xc <- Xc[idx]
    quo <- (Xc/Xm)
    if (uv_used == 'median'){
      d <- median(quo, na.rm = T)
    } else if (uv_used == 'mode'){
      den <- density(quo[!is.na(quo)], bw = width)
      d <- den$x[which.max(den$y)]
    }
    return(d)
  })
  Xn <- t(sapply(1:nrow(X), function(z){
      X[z,]/dilf[z]
  }))
  return(list(Xn = Xn, dilf = dilf))
}





