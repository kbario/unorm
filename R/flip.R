#' @title flipping spectra
#' @importFrom metabom8 get_idx

flip = function(X, ppm, i, sh=c(3, 3.1)){
  iid= get_idx(sh, ppm)
  out=apply(X, 1, function(x, idx=iid){
    if (sum(i[idx])<0) {
      i=i*-1
    }
    return(x)
  })
  t(out)
}
