#pqn alt

pqnam <- function(X, ppm, noi, use_ta = F, uv_used = 'mode', width){
  if (use_ta){
    cat('\033[0;34mPerforming Total Area Normalisation... ')
    X <- t(sapply(1:nrow(X), function(x){
      (X[x,])/(sum(X[x,]))
    }))
    cat('\033[1;32mDone.\n')
  }
  idx <- get_idx(c(0.25,9.5), ppm)
  cat('\033[0;34mCreating Reference Spectra... ')
  Xm <- apply(X, 2, median)
  Xm <- Xm[idx]
  cat('\033[1;32mDone.\n')
  cat('\033[0;34mCalculating Dilfs... ')
  dilf <- sapply(1:nrow(X), function(y){
    Xc <- X[y,idx]
    n <- noi[y]
    Xc[Xc<=n] = NA
    quo <- (Xc/Xm)
    if (uv_used == 'median'){
      d <- median(quo, na.rm = T)
    } else if (uv_used == 'mode'){
      den <- density(quo[!is.na(quo)], bw = width)
      d <- den$x[which.max(den$y)]
    }
    return(d)
  })
  cat('\033[1;32mDone.\n')
  cat('\033[0;34mNormalising X... ')
  Xn <- t(sapply(1:nrow(X), function(z){
    X[z,]/dilf[z]
  }))
  cat('\033[1;32mDone.\n')
  return(list(Xn = Xn, dilf = dilf))
}

