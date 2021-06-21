#' @title flipping spectra
#' @description NMR urine spectra that are produced from a small number of scans are often incorrectly orientated due to the water signal. This function looks at the creatinine signal (ppm = 3.05) which is in all urine spectra and see if the value is positive or negative. For negative valued spectra, they are mulitplied by -1 and flipped to the correct orientation.
#' @importFrom metabom8 get_idx

flip = function(X, ppm, sh=c(3, 3.1)){
  iid = get_idx(sh, ppm)
  out = t(apply(X, 1, function(x, idx=iid){
    if (sum(x[idx])<0) {
      x=x*-1
    }
    return(x)
  }))
}

