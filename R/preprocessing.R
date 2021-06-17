#' @title Standard 1D NMR Preprocessing
#' @details This function streamlines the preprocessing of NMR urine spectra by combining a range of functions that:
#' 1. orientate the spectra correctly,
#' 2. calibrated the spectra by a specific peak,
#' 3. calculate the line widths of the peaks and returns a warning with the spectra that exceed the specified threshold,
#' 4. remove the lower, upper, water and urea regions of the spectra,
#' 5. correct the baseline of the spectra using asymmetric least squares
#' 6. verify that the resulting X, ppm and meta objects match appropriately.
#' @param ppm An array of the chemical shift variables the same length as X.
#' @param X A matrix containing the NMR spectral data, the rows containing all values of a single experiment and the columns contain values per chemical shift variables.
#' @param meta The matrix of metadata pertaining to the X matrix. This is crucial for the TSP calibration and linewidth calculation.
#' @param flip A function that checks the orientation of spectra. This is a particularly important function for spectra derived from a small number of scans; the NMR orientates spectra based on what side of the y-axis has the most signal and water has a very large impact on this and without a large amount of metabolite signal on the correct side of the spectra (as seen in spectra derived from small scans) the spectra will be orientated the wrong way.
#' @param cali This is a metabom8 function that calibrates the spectra based on the TSP signal. Calibration is standardly performed by the NMR but this can be performed again for certainty.
#' @param calibrant This is the signal you wish to use as the calibrant. This signal must be compatible with the metabom8 function. For urine spectra, tsp is the main calibrant.
#' @param lineWid The maximum line width of peaks that should be filtered for. 1.0 is default and standard in NMR spectra.
#' @param baseline This function calls on the metabom8 function *bcor* to flatten and correct the baseline; baseline correction.
#' @param lowerCutoff The ppm value from where the lowest limit of the spectra will be trimmed/cutoff at.
#' @param upperCutoff The ppm value from where the highest limit of the spectra will be trimmed/cutoff at.
#' @param waterLower The lower ppm value from which the water region will be cut from.
#' @param waterUpper The upper ppm value from which the water region will be cut to.
#' @param ureaLower The lower ppm value from which the urea region will be cut from.
#' @param ureaUpper The upper ppm value from which the urea region will be cut to.
#' @importFrom metabom8 calibrate lw get_idx bcor
#' @export
#' @author \email{kylebario1@@gmail.com}
#' @family preproc
#' @examples
#' preProcessed <- preprocessing(X, ppm, meta, 0.5, 9.5, 4.7,4.85, 5.6,6)
#' preProcessed <- preprocessing(X, ppm, meta, flip = FALSE, cali = FALSE, baseline = TRUE, 0.5, 9.5, 4.7,4.85, 5.6,6)

preprocessing <- function(X, ppm, meta, flip = TRUE, cali = TRUE, calibrant = 'tsp', lineWid = 1.0, baseline = TRUE, lowerCutoff, upperCutoff, waterLower, waterUpper, ureaLower, ureaUpper){
  #relabel X and ppm
  X0 <- X
  ppm0 <- ppm

  #flip disorientated spectra
  if (flip){
    cat("flip was set to TRUE, flipping the spectra... ")
    Xf <- flip(X, sh=c(2.9,3))
    cat('Done\n')
  } else {
    warning("flip was set to FALSE, spectra have not been checked for orientation")
    Xf <- X
  }

  #calibrate spectra to tsp
  if (cali){
    cat("cali is set to TRUE, calibrating to ", calibrant,"... ")
    Xc <- calibrate(Xf, ppm, type = calibrant)
    cat('Done.\n')
  } else {
    cat('cali is set to FALSE...\n')
    cat('checking the metadata to verify calibration not required...\n')
    if (all(meta$p_SREF_mod==1)){
      Xc <- Xf
      cat('all spectra have already been calibrated.\n')
    } else {
      warning('not all spectra are calibrated, beginning calibration... ')
      Xc <- calibrate(Xf, ppm, type = 'tsp')
      print('Done.\n')
    }
  }

  # calculate the line width as quality control
  cat('checking linewidth of spectral peaks...\n')
  lwd <- lw(Xc,ppm,shift = c(-0.1,0.1))*meta$a_SFO1
  lwB <- lwd<lineWid
  if (all(lwB)){
    cat('all spectra have linewidths less than ', lineWid, '\n')
  } else {
    warning('the spectra ', which(lwB==FALSE), ' have linewidths greater than ',lineWid, ', check pproc[[1]] for further information')
  }
  DfX <- data.frame(lwd = lwd, lwB = lwB)
  DfX$names <- 1:(length(rownames(X)))
  DfX$orig_names <- rownames(X)

  #remove regions
  cat('removing the lower, water, urea and upper non-quantative regions of the spectra... ')
  idx_rm <- c(get_idx(range = c(min(ppm), lowerCutoff), ppm), get_idx(range = c(waterLower,waterUpper), ppm), get_idx(range = c(ureaLower, ureaUpper), ppm), get_idx(range = c(upperCutoff, max(ppm)), ppm))
  # remove these indexes
  Xr <- Xc[,-idx_rm]
  ppm <- ppm[-idx_rm]
  cat('Done.\n')

  # baseline correct the spectra
  if (baseline){
    cat('performing baseline correction... ')
    X <- bcor(Xr)
    cat('Done\n')
  }else{
    warning('baseline was set to FALSE. no baseline correction has been performed.')
    X <- Xr
  }


  #check meta and X have same rows
  cat('checking that the rows of X and meta match...\n')
  if (dim(meta)[1]==dim(X)[1]){
    cat('the number of spectra in X and meta match')
  } else {
    stop('the number of spectra in X and meta do not match')
  }
  #check X and ppm have same columns
  cat('checking that the length of ppm matches that of X...\n')
  if (dim(X)[2]==length(ppm)){
    cat('ppm and X columns match')
  } else {
    stop('X and ppm columns do not match')
  }
  pproc <- list(DfX, X0,ppm0, Xf, Xc, Xr, Xbl, X,ppm)
  return(pproc)
}
