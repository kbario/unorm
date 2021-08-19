#' Peak-Picking
#'
#' @param X The X array or matrix you want to pick the peak of
#' @param ppm The column matched ppm variable to your X
#' @param sh The chemical **sh**ift region you want to search
#' @param pm The amount you want to add to each side of the peak
#'
#' @return The ppm variables that match the peak you are searching for
#' @export
#'
#' @examples
shift_pickr <- function(X, ppm = NULL, sh = c(3,3.1), pm = 0.005){
  if (is.null(dim(X)) & is.null(length(X))){
    stop("X is neither a spectrum or a matrix, please provide an approproiate X.")
  }
  if (is.null(dim(X)) & !is.null(length(X))){
    if (is.null(ppm)){
      stop('Please provide a X-matched ppm variable')
    }
    p <- ppm
    if (length(sh)>2){
      i <- sh
    } else {
      i <- get_idx(sh, p)
    }
    m <- unname(which(X==max(X[i])))
    s <- get_idx(c(p[m]-pm,p[m]+pm), p)
    return(s)
  }
  if (!is.null(nrow(X))){
    if (is.null(ppm)){
      p <- as.numeric(colnames(X))
    } else {
      p <- ppm
    }
    if (length(sh)>2){
      i <- sh
    } else {
      i <- get_idx(sh, p)
    }
    s <- t(sapply(1:nrow(X), function(y){
      m <- unname(which(X[y,]==max(X[y,i])))
      s <- get_idx(c(p[m]-pm,p[m]+pm), p)
      return(s)
    }))
    return(s)
  }
}
