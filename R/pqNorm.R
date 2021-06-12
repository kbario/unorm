#' @title Probabilistic Quotient Normalisation
#' @param X a numerical matrix with rows being the spectra and columns being the chemical shift variables
#' @return a list of the normalised X matrix and a numerical array of the corresponding dilution factors
#' @details https://doi.org/10.1021/ac051632c
#' @author \email{kylebario1@@gmail.com}
#' @example
#' Xn <- pqn(X)
#' @export
#' @importFrom stats median density

pqNorm <- function(X, pcent = F){
  Xm <- apply(X, 2, median)
  Xm[Xm <= 0]=1
  dilf <- sapply(1:nrow(X), function(x){
    Xt <- X[x,]
    Xt[Xt<=0]=1
    quo <- (Xt/Xm)
    hgram <- hist(quo)
    vic <- c(min(hgram$breaks), seq(from = -0.005, to = 5.005, length.out = 502), max(hgram$breaks))
    h <- hist((log10(abs((Xt/Xm)))), breaks = vic, xlim= c(0.9,1.1))
    d <- h$mids[(which.max(h$counts[c(2:501)])+1)]
    return(d)
  })

  Xn <- sapply(1:nrow(X), function(y){
    if(dilf[y] == 1){
      X[y,]
    } else {
      X[y,]/dilf[y]
    }
  })
  if (pcent){
    p_cent <- sapply(1:nrow(X), function(x){
      hgram <- hist((X[x,]/Xm))
      vic <- c(min(hgram$breaks), seq(from = -0.005, to = 5.005, length.out = 502), max(hgram$breaks))
      h <- hist((X1[x,]/Xm), breaks = vic)
      p_cent <- ((sum(h$counts[c(2:501)]))/(sum(h$counts)))*100
      return(p_cent)
    })
    return(list(Xn, dilf, p_cent))
  }else{
    return(list(Xn, dilf))
  }
}



pqNorm(X1)

