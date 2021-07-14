#binning data by 0.01 ppm

binin <- function(X, ppm, bin_wit = 0.01){
  haf_bin <- 0.01/2
  low <- ppm_lim[1]-haf_bin
  upp <- ppm_lim[2]+haf_bin
  ppm_range <- upp-low
  n_bins <- ppm_range/bin_wit
  cuts <- seq(from = low, to = upp, by = bin_wit)
  idx <- sapply(1:n_bins, function(i){
    j <- i+1
    id <- get_idx(c(cuts[i],cuts[j]),ppm)
    return(id)
  })
  mids <- seq(from = ppm_lim[1], to = ppm_lim[2], by = bin_wit)
  Xb <- t(sapply(1:nrow(X),function(m){
    sapply(1:n_bins, function(n){
      ind <- unlist(idx[n])
      s <- sum(X[m,ind])
      return(s)
    })
  }))
  return(list(Xb = Xb, ppm_mids = mids))
}

npoints= 1030
iid = floor(length(ppm)/npoints)
iid = floor(length(ppm_new)/step)
ybin = rep(seq(iid), each = step)
bed <- approxfun(ppm, X[1,])
bed(ppm_new)
