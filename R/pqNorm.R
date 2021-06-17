#' @title Probabilistic Quotient Normalisation
#' @param X a numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @return A list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @details PQN is currently considered the gold standard NMR spectra normalisation method at this point in time. Multiple studies have tested its validity. PQN reliably normalises spectra and handles large amounts of noise. One limitation of the method is that it operates on the assumption that over 50% of the data must stay reasonably constant so hugely varying spectra may produces suboptimally normalised spectra. The methods paper can be found here: https://doi.org/10.1021/ac051632c
#' @author \email{kylebario1@@gmail.com}
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
