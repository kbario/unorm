#' @title Standard 1D NMR Preprocessing
#' @description [Preproc()] is a function aimed at streamlining the preprocessing stage of analysing NMR spectra. It harnesses the power of the package `metabom8` to quickly and easily process spectra.
#' @details [Preproc()] is powered by the `metabom8` package and is simply a tool to harness and streamline functions in `metabom8` when preprocessing NMR spectra, not replace them.
#' # The Pipeline
#' This function streamlines the preprocessing of NMR urine spectra by combining a range of functions. It:
#' 1. Orientates the spectra correctly,
#' 2. Calibrates the spectra by a specific peak,
#' 3. Calculates the line widths of the peaks and returns a warning with the spectra that exceed the specified threshold,
#' 4. Removes the lower, upper, water and urea regions of the spectra,
#' 5. Corrects the baseline of the spectra using asymmetric least squares
#' 6. Verifies that the resulting X, ppm and meta objects match appropriately.
#' @param ppm An array of the chemical shift variables, column matched to X.
#' @param X A matrix containing the non-preprocessed NMR spectral data. The rows should containing all values of a single experiment, and the columns, the values of the chemical shift variables.
#' @param meta The matrix of metadata pertaining to the X matrix. This is crucial for the TSP calibration and line width calculation.
#' @param flip Default is set to `FALSE`. **NOT CRUCIAL FOR STANDARDLY ACQUIRED SPECTRA** This function checks the orientation of spectra. This is a particularly important function for spectra derived from a **small number of scans**; the NMR orientates spectra based on what side of the y-axis has the most signal and water has a very large impact on this and without a large amount of metabolite signal on the correct side of the spectra (as seen in spectra derived from small scans) the spectra will be orientated the wrong way.
#' @param cali Default is set to `FALSE`. This calls on the `metabom8` function [metabom8::calibrate()] which ensures peaks are aligned  based on the TSP signal. Calibration is standardly performed by the NMR but this can be performed again for certainty.
#' @param calibrant This is the signal you wish to use as the calibrant. This signal must be compatible with [metabom8::calibrate()]. For urine spectra, *tsp* is the main calibrant.
#' @param lineWid This argument provides import for [metabom8::lw()], a function that calculates the line width of peaks. Spectra with peaks that have small line widths have sharper and more precise results which is more desirable. A maximum cutoff of 1 ensures spectra contain robust results. Consider omitting spectra with line widths over 1.
#' @param baseline Default is set to `TRUE`. This argument calls on [metabom8::bcor()], a baseline correcting function to smooth the spectral baselines and remove the influence of broad peaks.
#' @param lowCut A single floating point number defining the ppm value that the lower limit of the spectra are trimmed to.
#' @param uppCut A single floating point number defining the ppm value that the upper limit of the spectra are trimmed to.
#' @param watCut The lower and upper ppm values concatenated, from which the water region will be trimmed and omitted. Water regions provide no important information and should be removed prior to data analysis. Default is set to `c(4.5,5)`
#' @param ureCut The lower and upper ppm values concatenated, from which the urea region will be trimmed and omitted. Urea regions also provide no important information and should be removed prior to data analysis. Default is set to `c(5.6,6)`
#' @param noi_sh The shift of the noise region used to calculate estimation of noise
#' @importFrom metabom8 calibrate lw get_idx
#' @return This function returns a list with:
#' 1. The processed X matrix in the first element,
#' 2. The processed ppm array in the second element, and
#' 3. The line width results in a data frame in the third element.
#' * Following the example below will extract the results quickly and easily.
#' @export
#' @author \email{kylebario1@@gmail.com}
#' @family {preproc}
#' @examples
#' path = system.file('extdata', package = 'unorm')

Preproc <- function(X, ppm, meta, baseline = T, flip = F, cali = F, calibrant = 'tsp', lineWid = 1.0, lowCut = 0.25, watCut = c(4.5,5), ureCut = c(5.6,6), uppCut = 9.5, noi_sh = c(9.5, 11)){
  #relabel X and ppm
  X_OG <- X
  ppm_OG <- ppm

  #flip disorientated spectra
  if (flip){
    cat('\033[0;34mFlipping the spectra... \033[0m')
    Xf <- flip(X, ppm)
    cat('\033[1;32mDone.\n\033[0m')
  } else {
    cat("\033[1;33mSpectra haven't been checked for orientation.\n\033[0m")
    Xf <- X
  }

  #calibrate spectra to tsp
  if (cali){
      cat('\033[0;34mCalibrating to ', calibrant,"... ", sep = '')
      Xc <- calibrate(Xf, ppm, type = calibrant)
      cat('\033[1;32mDone.\n\033[0m')
    } else {
    cat('\033[0;34mChecking that calibration isn\'t required... \033[0m')
    if (all(meta$p_SREF_mod==1)){
      Xc <- Xf
      cat('\033[1;32mThey\'re all good.\n\033[0m')
    } else {
      cat('\033[1;33mCalibration is required\n\033[0m')
      cat('\033[0;34mCalibrating to', calibrant,"... ")
      Xc <- calibrate(Xf, ppm, type = calibrant)
      cat('\033[1;32mDone.\n\033[0m')
    }
  }

  # calculate the line width as quality control
  cat('\033[0;34mChecking line width of spectra... \033[0m')
  lwd <- lw(Xc,ppm,shift = c(-0.1,0.1), sf = meta$a_SFO1)
  lwB <- lwd<lineWid
  if (all(lwB)){
    cat('\033[1;32mAll spectra have linewidths under', lineWid, '\n\033[0m')
  } else {
    cat('\033[1;31m', length(which(lwB==FALSE)), 'spectra have line-widths over',lineWid,'.\n')
    cat('\033[1;33mCheck Df_ppro for more information.\n\033[0m')
  }
  DfX <- data.frame(lwd = lwd, lwB = lwB)
  DfX$names <- 1:(length(rownames(X)))
  DfX$orig_names <- rownames(X)

  #remove regions
  cat('\033[0;34mRemoving non-quantative regions... \033[0m')
  idx_rm <- c(get_idx(range = c(min(ppm), lowCut), ppm), get_idx(range = watCut, ppm), get_idx(range = ureCut, ppm))
  # remove these indexes
  Xr <- Xc[,-idx_rm]
  ppm <- ppm[-idx_rm]
  cat('\033[1;32mDone.\n\033[0m')

  # baseline correct the spectra
  if (baseline){
    cat('\033[0;34mPerforming baseline correction... \033[0m')
    X <- bl(Xr)
    cat('\033[1;32mDone.\n\033[0m')
  }else{
    cat('\033[1;33mNo baseline performed.\n\033[0m')
    X <- Xr
  }

  #calculating the noise
  cat('\033[0;34mCalculating Noise Estimations... \033[0m')
  noi <- noise(X, ppm, shift = c(noi_sh[1], noi_sh[2]), sd_mult = 5)
  cat('\033[1;32mDone.\n\033[0m')

  X <- X[,-get_idx(range = c(uppCut, max(ppm)), ppm)]
  ppm <- ppm[-get_idx(range = c(uppCut, max(ppm)), ppm)]

  #check meta and X have same rows
  cat('\033[0;34mChecking that X and meta rows match... \033[0m')
  if (dim(meta)[1]==dim(X)[1]){
    cat('\033[1;32mDone.\n\033[0m')
  } else {
    stop('The number of spectra in X and meta do not match.\n')
  }
  #check X and ppm have same columns
  cat('\033[0;34mChecking that ppm length and X columns match... \033[0m')
  if (dim(X)[2]==length(ppm)){
    cat('\033[1;32mDone.\n\033[0m')
  } else {
    stop('X columns and ppm do not match.\n')
  }
  assign("X", X, envir = .GlobalEnv)
  assign("ppm", ppm, envir = .GlobalEnv)
  assign("Df_ppro", DfX, envir = .GlobalEnv)
  assign("noi", noi, envir = .GlobalEnv)
  assign("X_OG", X_OG, envir = .GlobalEnv)
  assign("ppm_OG", ppm_OG, envir = .GlobalEnv)
}
