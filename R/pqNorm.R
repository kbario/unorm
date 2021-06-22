#' @title Probabilistic Quotient Normalisation
#' @description PQN is currently the gold standard method used to normalise NMR spectra.
#' @details
#' ### How It Works:
#' PQN works by normalising experimental spectra in the provided X matrix in relation to a reference spectrum.
#' 1. [pqNorm()] creates the reference automatically by calculating the median spectrum of X as outlined in the initial methods paper (see 'See also')
#' 2. [pqNorm()] derives a quotient for each ppm value within the limits of shift in a experimental spectrum by dividing it's intensities with that of the reference spectrum's.
#' 3. The most frequently occurring (median) quotient is calculated and is said to be the dilution coefficient (`dilf`) of that spectrum.
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
#' @return The output of this function is a list containing:
#' 1. The normalised version of X in the first element and
#' 2. A numerical array of the corresponding dilution factors calculated by the function.
#' * Following the example below will extract the results quickly and easily.
#' @seealso The methods paper first describing PQN can be found here: \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' Xpq <- pqNorm(X)
#' Xpqn <- Xpq[[1]])
#' pqnDilf <- Xpq[[2]]
#' @importFrom metabom8 bcor get_idx
#' @importFrom graphics hist
#' @export

pqNorm <- function(X, ppm, shift = c(0.5,8.5)){
  Xta <- t(sapply(1:nrow(X), function(x){
    (X[x,])/(sum(X[x,]))
  }))
  idx <- get_idx(shift, ppm)
  Xc <- Xta[, idx]
  Xm <- apply(Xc, 2, median)
  dilf <- sapply(1:nrow(X), function(y){
    Xt <- Xc[y,]
    quo <- (Xt/Xm)
    d <- median(quo)
    return(d)
  })
  Xn <- t(sapply(1:nrow(X), function(z){
      X[z,]/dilf[z]
  }))
  return(list(Xn, dilf))
}





