<<<<<<< HEAD
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
=======
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
#' @param linWid The maximum line width of peaks that should be filtered for. 1.0 is default and standard in NMR spectra.
#' @param baseline This function calls on the metabom8 function *bline* to flatten and correct the baseline; baseline correction.
#' @param lowerCutoff The ppm value from where the lowest limit of the spectra will be trimmed/cutoff at.
#' @param upperCutoff The ppm value from where the highest limit of the spectra will be trimmed/cutoff at.
#' @param waterLower The lower ppm value from which the water region will be cut from.
#' @param waterUpper The upper ppm value from which the water region will be cut to.
#' @param ureaLower The lower ppm value from which the urea region will be cut from.
#' @param ureaUpper The upper ppm value from which the urea region will be cut to.
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
#' @importFrom metabom8 calibrate lw get_idx bline
#' @example
#'     preprocessed <- preprocessing(X, ppm, meta, flip = F, cali = F, 0.5, 9.5, 4.5, 5, 5.2, 6)
#' @export
#' @family preproc
#' @examples
#' preProcessed <- preprocessing(X, ppm, meta, 0.5, 9.5, 4.7,4.85, 5.6,6)
#' preProcessed <- preprocessing(X, ppm, meta, flip = FALSE, cali = FALSE, baseline = TRUE, 0.5, 9.5, 4.7,4.85, 5.6,6)

<<<<<<< HEAD
preprocessing <- function(X, ppm, meta, flip = TRUE, cali = TRUE, caliType = 'tsp', baseline = TRUE, lowerCutoff, upperCutoff, waterLower, waterUpper, ureaLower, ureaUpper){
=======
preprocessing <- function(X, ppm, meta, flip = TRUE, cali = TRUE, calibrant = 'tsp', linWid = 1.0, baseline = TRUE, lowerCutoff, upperCutoff, waterLower, waterUpper, ureaLower, ureaUpper){
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  #relabel X and ppm
  X0 <- X
  ppm0 <- ppm

  #flip disorientated spectra
  if (flip){
<<<<<<< HEAD
    cat('spectra are being flipped... ')
    Xf <- flip(X, sh=c(2.9,3))
    cat('Done\n')
  } else {
    cat('flip = FALSE, spectra have not been checked for orientation\n')
    Xf <- X
=======
    cat("flipping the spectra... ")
    Xf <- flip(X, sh=c(2.9,3))
    cat('Done\n')
  } else {
    warning("flip was set to FALSE, spectra have not been checked for orientation")
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  }

  #calibrate spectra to tsp
  if (cali){
<<<<<<< HEAD
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
=======
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
  lwB <- lwd<linWid
  if (all(lwB)){
    cat('all spectra have linewidths less than ', linWid, '\n')
  } else {
    warning('the spectra ', which(lwB==FALSE), ' have linewidths greater than ',linWid, ', check pproc[[1]] for further information')
  }
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  DfX <- data.frame(lwd = lwd, lwB = lwB)
  DfX$names <- 1:(length(rownames(X)))
  DfX$orig_names <- rownames(X)
  if (lwB){
    cat('all spectra have a linewidth less than 1.0\n')
  } else{
    warning('the spectra ', which(DfX$lwB==FALSE), ' have linewidths over 1.0. check pproc[[1]] for further information')
  }

  #remove regions
<<<<<<< HEAD
  cat('removing the lower, water, urea and upper non-quantative regions of the spectra... ')
=======
  cat('removing lower, water, urea and upper regions... ')
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  idx_rm <- c(get_idx(range = c(min(ppm), lowerCutoff), ppm), get_idx(range = c(waterLower,waterUpper), ppm), get_idx(range = c(ureaLower, ureaUpper), ppm), get_idx(range = c(upperCutoff, max(ppm)), ppm))
  # remove these indexes
  Xr <- Xc[,-idx_rm]
  ppm <- ppm[-idx_rm]
<<<<<<< HEAD
  cat('Done\n')
=======
  cat('Done.\n')
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819

  # baseline correct the spectra
  if (baseline){
    cat('performing baseline correction... ')
    Xbl <- bline(Xr)
<<<<<<< HEAD
    print('Done\n')
    X <- Xbl
  }else{
    X <- Xr
    warning('baseline = FALSE, no baseline was performed')
=======
    cat('Done.')
  } else {
    warning('baseline was set to FALSE. no baseline correction has been performed.')
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  }


  #check meta and X have same rows
  cat('checking that the rows of X and meta match...\n')
  if (dim(meta)[1]==dim(X)[1]){
    cat('the number of spectra in X and meta match')
  } else {
<<<<<<< HEAD
    error('the number of spectra in X and meta do not match')
=======
    stop('the number of spectra in X and meta do not match')
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  }
  #check X and ppm have same columns
  cat('checking that the length of ppm matches that of X...\n')
  if (dim(X)[2]==length(ppm)){
    cat('ppm and X columns match')
  } else {
<<<<<<< HEAD
    error('X and ppm columns do not match')
=======
    stop('X and ppm columns do not match')
>>>>>>> 111289ecd790675ab4e3930d075cb8825f040819
  }
  pproc <- list(DfX, X0,ppm0, Xf, Xc, Xr, Xbl, X,ppm)
  return(pproc)
}
