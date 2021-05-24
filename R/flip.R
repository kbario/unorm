#' @title flipping spectra
#' @importFrom metabom8 get_idx

flip = function(X, ppm, x, sh=c(3, 3.1)){
  iid= get_idx(sh, ppm)
  out=apply(X, 1, function(x, idx=iid){
    if (sum(x[idx])<0) {
      x=x*-1
    }
    return(x)
  })
  t(out)
}
