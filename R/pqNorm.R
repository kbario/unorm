#' @title Probabilistic Quotient Normalisation (PQN)
#' @description [pqNorm()] is `r packageName()`'s function to normalise spectra using PQN.
#' @param X The numerical matrix containing the NMR data. The rows contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @return The output of this function is a list containing the normalised version of X and a numerical array of the corresponding dilution factors.
#' @details ## How It Works:
#' * PQN works by normalising experimental spectra in the provided X matrix in relation to a reference spectrum.
#' * [pqNorm()] creates the reference automatically by calculating the median spectrum of X as outlined in the initial methods paper (see 'See also')
#' * [pqNorm()] derives a quotient for each ppm value in a experimental spectrum by dividing it's intensities with that of the reference spectrum's.
#' * The most frequently occurring quotient is calculated and is said to be the dilution coefficient (`dilf`) of that spectrum.
#' * The sample is then scaled with this dilf and will be comparable with all other spectra normalised with the reference spectrum.
#' ## Background Information:
#' * PQN is currently the gold standard for normalising NMR spectra.
#' * Multiple studies have tested its validity.
#' * PQN reliably normalises spectra and handles large amounts of noise.
#' * PQN's biggest limitation is that it operates on the assumption that the majority of the spectra will stay constant so hugely varying spectra may be suboptimally normalised.
#' @seealso The methods paper first describing PQN can be found here: \url{https://doi.org/10.1021/ac051632c}
#' @author \email{kylebario1@@gmail.com}
#' @family {Reference-Based}
#' @examples
#' Xn <- pqNorm(X)
#' pqnDilf <- Xn[[2]]
#' Xn <- t(Xn[[1]])
#' @importFrom metabom8 bcor
#' @importFrom graphics hist
#' @export

pqNorm <- function(X){
  Xm <- apply(X, 2, median)
  Xm <- bcor(Xm)
  Xm[Xm <= 0]=0
  dilf <- sapply(1:nrow(X), function(x){
    Xt <- X[x,]
    Xt[Xt<=0]=0
    quo <- (Xt/Xm)
    hgram <- hist(quo, plot = F)
    vic <- c(min(hgram$breaks), seq(from = 0.005, to = 5.005, length.out = 501), max(hgram$breaks))
    h <- hist(Xt/Xm, breaks = vic, xlim = c(0,5))
    d <- h$mids[(which.max(h$counts[c(2:501)])+2)]
    d
    return(d)
  })
  Xn <- sapply(1:nrow(X), function(y){
      X[y,]/dilf[y]
  })
  return(list(Xn, dilf))
}
