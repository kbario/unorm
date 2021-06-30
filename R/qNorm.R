qNorm <- function(X){
  Xs <- t(apply(X, 1, sort))
  Xr <- t(apply(X, 1, rank))
  Xm <- apply(Xs, 2, mean)
  Xq <- t(sapply(1:nrow(X), function(x){
    Xm[Xr[x,]]
  }))
  Xsum <- t(apply(X, 1, sum))
  Xqsum <- t(apply(Xq, 1, sum))
  dilf <- sapply(1:nrow(X), function(h){
    Xsum[h]/Xqsum[h]})
  Xn <- t(sapply(1:nrow(X), function(x){
    X[x,]/dilf[x]
  }))
  return(list(Xn = Xn, dilf = dilf))
}

Xq <- qNorm(X)

binning <- function(X, ppm, ppm_range = c(0.5, 9.5), width = 0.1, hop = 'ppm', npoint){
  s <- ppm_range[1]
  nbins <- ((ppm_range[2])-(ppm_range[1]))/width
  sq <- seq(from = ppm_range[1], to = ppm_range[2], by = width)
  idxs <- sapply(1:nbins, function(x){
    y = x+1
    i <- get_idx(c(sq[x], sq[y]), ppm)
    return(i)
  })
  Xb <- t(sapply(1:nrow(X), function(x){
    sapply(1:nbins, function(y){
      z <- idxs[[y]]
      sum(X[x,z])
    })
  }))
  return(Xb)
}
Xb <- binning(X, ppm, ppm_range = c(0.5,9.5), width = 0.1)


Xm <- apply(Xs, 2, mean)





aa <- c(1:9)
bb <- c(9:1)
aa[bb]
aa <- Xo[1,]
Xm[aa]
Xo[1,Xm]
















