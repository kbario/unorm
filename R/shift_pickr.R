#' Peak-Picking
#'
#' @param X The X array or matrix you want to pick the peak of
#' @param ppm The column matched ppm variable to your X
#' @param sh The chemical **sh**ift region you want to search
#' @param pm The amount you want to add to each side of the peak
#'
#' @return The ppm variables that match the peak you are searching for
#' @export
#' @family {estimation}
#'
#' @examples
#' data(X)
#' idx <- shift_pickr(X[1,], ppm, sh = c(4,4.1), pm = 0.005)
#' idx1 <- shift_pickr(X[2,], ppm, sh = c(4,4.1), pm = 0.01)
#' plot(ppm[idx], X[1,idx], main = "Creatinine Peak 4.05", type = 'l', xlim = c(4.05, 4.07))
#' points(ppm[idx1], X[2,idx1], type = 'l', col = 'red')
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
