#' @title Standard Preprocessing of 1D spectra
#' @param X a numerical matrix with rows being the spectra and the columns being the chemical shift variables
#' @param ppm a numerical array holding the chemical shift values of the X matrix
#' @param meta a matrix containing the metadata of X
#' @param flip specifies if the flip function will correctly orientate the spectra. highly recommended for spectra with a low number of scans. default is TRUE
#' @param cali specifies if the metabom8 function 'calibrate' is performed. this is not necessary if this has been done at an NMR level. default is TRUE
#' @param caliType specifies what signal will be used to calibrate the spectra. default is 'tsp'
#' @param baseline specifies if the metabom8 function 'baseline' will be performed. highly recommended this is performed. default = TRUE
#' @param lowerCutoff requires a numerical value that defines the ppm value from where the spectra will be cutoff at the lower end
#' @param upperCutoff requires a numerical value that defines the ppm value from where the spectra will be cutoff at the upper end
#' @param waterLower a numerical value of the ppm defining the lower cutoff of the water
#' @param waterUpper a numerical value of the ppm defining the upper cutoff of the water
#' @param ureaLower a numerical value of the ppm defining the lower cutoff of the urea
#' @param ureaUpper a numerical value of the ppm defining the upper cutoff of the urea
#' @author \email{kylebario1@@gmail.com}
#' @details this function will take a matrix containing spectra in the rows and chemical shift values in the columns, the a numerical array of chemical shifts, and the metadata of the X matrix, and ensure that it is correctly orientated, calibrated, baseline corrected, and the noise regions removed, while also providing information about the linewidths of the peaks. errors will appear if X, ppm and meta do not match at the end of the function. contact author if this occurs.
#' @importFrom metabom8 calibrate lw get_idx bline
#' @export
#' @family preproc
#' @examples
#' preProcessed <- preprocessing(X, ppm, meta, 0.5, 9.5, 4.7,4.85, 5.6,6)
#' preProcessed <- preprocessing(X, ppm, meta, flip = FALSE, cali = FALSE, baseline = TRUE, 0.5, 9.5, 4.7,4.85, 5.6,6)

preprocessing <- function(X, ppm, meta, flip = TRUE, cali = TRUE, caliType = 'tsp', baseline = TRUE, lowerCutoff, upperCutoff, waterLower, waterUpper, ureaLower, ureaUpper){
  #relabel X and ppm
  X0 <- X
  ppm0 <- ppm

  #flip disorientated spectra
  if (flip){
    cat('spectra are being flipped... ')
    Xf <- flip(X, sh=c(2.9,3))
    cat('Done\n')
  } else {
    cat('flip = FALSE, spectra have not been checked for orientation\n')
    Xf <- X
  }

  #calibrate spectra to tsp
  if (cali){
    cat('checking if spectra require calibration... ')
    if (all(meta$p_SREF_mod==1)){
      Xc <- Xf
      cat('Done, spectra did not require calibration\n')
    } else {
      Xc <- calibrate(Xf, ppm, type = caliType)
      print('Done, spectra required calibration\n')
    }
  } else {
    cat('cali = FALSE, spectra were not calibrated\n')
  }

  # calculate the line width as quality control
  cat('checking linewidth of spectra... ')
  lwd <- lw(Xc,ppm,shift = c(-0.1,0.1))*meta$a_SFO1
  table(lwd<1.0)
  lwB <- lwd<1.0
  DfX <- data.frame(lwd = lwd, lwB = lwB)
  DfX$names <- 1:(length(rownames(X)))
  DfX$orig_names <- rownames(X)
  if (lwB){
    cat('all spectra have a linewidth less than 1.0\n')
  } else{
    warning('the spectra ', which(DfX$lwB==FALSE), ' have linewidths over 1.0. check pproc[[1]] for further information')
  }

  #remove regions
  cat('removing the lower, water, urea and upper non-quantative regions of the spectra... ')
  idx_rm <- c(get_idx(range = c(min(ppm), lowerCutoff), ppm), get_idx(range = c(waterLower,waterUpper), ppm), get_idx(range = c(ureaLower, ureaUpper), ppm), get_idx(range = c(upperCutoff, max(ppm)), ppm))
  # remove these indexes
  Xr <- Xc[,-idx_rm]
  ppm <- ppm[-idx_rm]
  cat('Done\n')

  # baseline correct the spectra
  if (baseline){
    cat('performing baseline correction... ')
    Xbl <- bline(Xr)
    print('Done\n')
    X <- Xbl
  }else{
    X <- Xr
    warning('baseline = FALSE, no baseline was performed')
  }


  #check meta and X have same rows
  if (dim(meta)[1]==dim(X)[1]){
    print('the number of spectra in X and meta match')
  } else {
    error('the number of spectra in X and meta do not match')
  }
  #check X and ppm have same columns
  if (dim(X)[2]==length(ppm)){
    print('X and ppm columns match')
  } else {
    error('X and ppm columns do not match')
  }
  pproc <- list(DfX, X0,ppm0, Xf, Xc, Xr, Xbl, X,ppm)
  return(pproc)
}
