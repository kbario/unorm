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
#' Xn <- pqNorm(X)
#' Xpqn <- Xn[[1]]
#' pqnDilf <- Xn[[2]]
#' @importFrom metabom8 bcor get_idx
#' @importFrom graphics hist
#' @importFrom stats sd
#' @export

pqNorm <- function(X, ppm, shift = c(0.5,9.5), noise = c(9.5,11), use_ta = F, dilf_calc_method = 'combination', precision = 0.1){
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
    int_rm <- 10*sd(Xc[get_idx(noise, ppm)])+mean(Xc[get_idx(noise, ppm)], trim = .05)
    Xc[Xc<=int_rm]=NA
    Xc <- Xc[idx]
    quo <- (Xc/Xm)
    if (dilf_calc_method == 'median'){
      d <- median(quo, na.rm = T)
    } else if (dilf_calc_method == 'mode'){
      ends <- precision/2
      if (min(quo, na.rm = T)<ends){
        h <- hist(quo, breaks = c(0, seq(from = ends, to = (ceiling(max(quo, na.rm = T))+ends), by = precision)),plot = F)#, xlim = c(0,5))
        d <- h$mids[which.max(h$counts)]
      } else {
        h <- hist(quo, breaks = seq(from = ends, to = (ceiling(max(quo, na.rm = T))+ends), by = precision), plot = F)#, xlim = c(0,5))
        d <- h$mids[which.max(h$counts)]
      }
    } else if (dilf_calc_method == 'combination'){
      me <- median(quo, na.rm = T)
      if (min(quo, na.rm = T)<ends){
        h <- hist(quo, breaks = c(0, seq(from = ends, to = (ceiling(max(quo, na.rm = T))+ends), by = precision)),plot = F)#, xlim = c(0,5))
        mo <- h$mids[which.max(h$counts)]
      } else {
        h <- hist(quo, breaks = seq(from = ends, to = (ceiling(max(quo, na.rm = T))+ends), by = precision), plot = F)#, xlim = c(0,5))
        mo <- h$mids[which.max(h$counts)]
      }
      d <- (me+mo)/2
    }
    return(d)
  })
  Xn <- t(sapply(1:nrow(X), function(z){
      X[z,]/dilf[z]
  }))
  return(list(Xpqn = Xn, pqnDilf = dilf))
}





