#' @title flipping spectra
#' @description Spectra produced from experiments with small scans can often appear upside-down. This function corrects this.
#' @details NMR urine spectra that are produced from a small number of scans are often incorrectly orientated due to the water signal. This function looks at the creatinine signal (ppm = 3.05) which is in all urine spectra to see if the value is positive or negative. For negative valued spectra, they are mulitplied by -1 and flipped to the correct orientation.
#' @family {preproc}
#' @param X The numerical matrix containing the NMR data you wish to check for orientation. The rows must contain information of one whole spectrum and the columns contain the specific chemical shift variables.
#' @param ppm An array of chemical shift variables. ppm should be column matched to the X matrix you are normalising.
#' @param shift The concatenated ppm values that define the lower and upper bounds of the creatinine signal. Default is c(3,3.1)
#' @return Returns the original X matrix but with all values with the correct sign.
#' @author \email{kylebario1@@gmail.com}
#' @family {Data Minipulation}
#' @examples
#' path = system.file('extdata', package = 'unorm')
#' Xf <- flip(X, ppm)
#' @importFrom metabom8 get_idx
#' @export

flip = function(X, ppm, shift = c(3, 3.1)){
  iid = get_idx(shift, ppm)
  out = t(apply(X, 1, function(x, idx=iid){
    if (sum(x[idx])<0) {
      x=x*-1
    }
    return(x)
  }))
}

