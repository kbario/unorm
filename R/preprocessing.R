#' @title standard preprocessing 1D
#' @param ppm an array of the chemical shift variables
#' @importFrom metabom8 calibrate lw get_idx bline
#' @export
#' @family preproc

preprocessing <- function(X, ppm, meta, flip = TRUE, cali = TRUE, baseline = TRUE, lowerCutoff, upperCutoff, waterLower, waterUpper, ureaLower, ureaUpper){
  #relabel X and ppm
  X0 <- X
  ppm0 <- ppm

  #flip disorientated spectra
  if (flip){
    Xf <- flip(X, sh=c(2.9,3))
    print('spectra the right way up')
  }



  #calibrate spectra to tsp
  if (all(meta$p_SREF_mod==1)){
    Xc <- Xf
    print('no tsp calibration required')
  } else {
    Xc <- calibrate(Xf, ppm, type = 'tsp')
    print('tsp calibration performed')
  }

  # calculate the line width as quality control
  lwd <- lw(Xc,ppm,shift = c(-0.1,0.1))*meta$a_SFO1
  table(lwd<1.0)
  print(table(lwd<1.0))
  lwB <- lwd<1.0
  DfX <- data.frame(lwd = lwd, lwB = lwB)
  DfX$n <- 1:(length(rownames(X)))
  DfX$orig_n <- rownames(X)

  #remove regions
  idx_rm <- c(get_idx(range = c(min(ppm), lowerCutoff), ppm), get_idx(range = c(waterLower,waterUpper), ppm), get_idx(range = c(ureaLower, ureaUpper), ppm), get_idx(range = c(upperCutoff, max(ppm)), ppm))
  # remove these indexes
  Xr <- Xc[,-idx_rm]
  ppm <- ppm[-idx_rm]
  print('non-quantitative regions removed')

  # baseline correct the spectra
  Xbl <- bline(Xr)
  print('baseline correction performed')

  # rename X
  X <- Xbl

  #check meta and X have same rows
  if (dim(meta)[1]==dim(X)[1]){
    print('the number of spectra in X and meta match')
  } else {
    print('the number of spectra in X and meta do not match')
  }
  #check X and ppm have same columns
  if (dim(X)[2]==length(ppm)){
    print('X and ppm columns match')
  } else {
    print('X and ppm columns do not match')
  }
  pproc <- list(DfX, X0,ppm0, Xf, Xc, Xr, Xbl, X,ppm)
  return(pproc)
}
