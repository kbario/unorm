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
#' @param noi Takes an array that is row matched to the X matrix you are normalising with the values equaling the maximum noise estimation for each spectra respectively.
#' @param use_ta Requires a boolean `TRUE` or `FALSE` if total area normalisation should be performed on the spectra before PQN is.
#' @param uv_used PQN utilises finding the median or the mode, which are both *U*ni*v*ariate methods. Recognises either the string 'median' or 'mode' to instruct which method to use. Default = 'mode'
#' @param calc_region The lower and upper bounds of the spectrum that will be used to calculate the dilution coeficient
#' @param bin_width The width of the bin when the spectra are binned
#' @return The output of this function is a list containing:
#' 1. The normalised version of X in the first element and
#' 2. A numerical array of the corresponding dilution factors calculated by the function.
#' * Following the example below will extract the results quickly and easily.
#' @seealso The methods paper first describing PQN can be found here: \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @examples
#' data(X, noi)
#' pqNorm(X, noi)
#' cat(dilf_pqn)
#' @importFrom metabom8 get_idx
#' @importFrom stats sd density median
#' @export

pqNorm <- function(X, noi, use_ta = F, uv_used = 'mode', calc_region = c(0.5,9.5), bin_width = 0.01){
  if (is.null(use_ta)){
    stop("use_ta was left blank. Please specify if you would like to normalise the spectra with Total Area prior to PQN with T or F.")
  }
  if (is.null(uv_used)){
    stop("uv_used was left blank. Please specify which univariate variable you want to use to calculate the dilf. 'mode' and 'median' accepted.")
  }
  if (!uv_used == 'mode' & !uv_used == 'median'){
    stop("uv_used not specified correctly. Only 'mode' and 'median' are accepted.")
  }
  if (!names(noi)[1]==rownames(X)[1] | !all(!is.na(match(names(noi), rownames(X)))) | !all(diff(match(names(noi), rownames(X)))==1)){
    if (length(noi) < nrow(X)){
      stop("noi does not contain as many values as X has spectra. To continue use uNorm::noise() to calculate noise estimations for your X.")
    }
    cat("\033[1;33mThe provided X and noi arguements do not match. Attempting to match them now... \033[0m")
    ma <- match(rownames(X),names(noi))
    noi_matched <- noi[ma]
    if (all(names(noi_matched)==rownames(X))){
      assign('noi_matched', noi_matched, envir = .GlobalEnv)
      cat('\033[1;32mDone.\n\033[0;34mFind noi_matched in the global environment.\033[0m')
    }
    if (!all(names(noi_matched)==rownames(X))){
      cat('\n\033[0;31mX and noi could not be matched. Please match them before applying PQN.\n\033[0m')
    }
  }
  if (use_ta){
    cat('\033[0;34mPerforming Total Area Normalisation... \033[0m')
    X <- t(sapply(1:nrow(X), function(x){
      (X[x,])/(sum(X[x,]))
    }))
    cat('\033[1;32mDone.\n\033[0m')
  }
  p <- as.numeric(colnames(X))
  cat('\033[0;34mPreparing Spectra and Reference...\n\033[0m')
  idx <- get_idx(c(calc_region[1],calc_region[2]), p)
  cat('\033[0;34mSelecting ppm and Removing Noise... \033[0m')
  Xs <- t(sapply(1:nrow(X), function(i){
    n <- noi[i]
    Xs <- X[i,idx]
    Xs[Xs<n]=0
    return(Xs)
  }))
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mBinning... \033[0m')
  Xsb <- binin(Xs, p[idx], width = bin_width)
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mCalculating Reference Spectrum... \033[0m')
  Xm <- apply(Xsb, 2, median)
  cat('\033[1;32mDone.\n')
  cat('\033[0;34mCalculating Dilfs... \033[0m')
  q <- t(sapply(1:nrow(Xsb), function(j){
      Xsb[j,]/Xm
  }))
  if (uv_used == "mode"){
    cat('\033[0;34mUsing the Mode... \033[0m')
    dilf <- sapply(1:nrow(q), function(y){
      i <- q[y,]
      d <- suppressWarnings(log10(i)[!is.nan(i) & !is.infinite(i) & !is.na(i)])
      den <- suppressWarnings(density(d[!is.nan(d) & !is.infinite(d) & !is.na(d)]))
      dilf <- 10^(den$x[which.max(den$y)])
      return(dilf)
    })
  }
  if (uv_used == "median"){
    cat('\033[0;34mUsing the Median... \033[0m')
    dilf <- sapply(1:nrow(q), function(z){
      i <- q[z,]
      dilf <- median(i, na.rm = T)
      return(dilf)
    })#

  }
  cat('\033[1;32mDone.\n\033[0m')
  cat('\033[0;34mNormalising X... \033[0m')
  Xn <- t(sapply(1:nrow(X), function(a){
    X[a,]/dilf[a]
  }))
  rownames(Xn) <- rownames(X)
  cat('\033[1;32mDone.\n\033[0m')
  assign("X_pqn", Xn, envir = .GlobalEnv)
  assign("dilf_pqn", dilf, envir = .GlobalEnv)
  assign('X_pqn_median', Xm, envir = .GlobalEnv)
  assign('X_pqn_binned', Xsb, envir = .GlobalEnv)
}
