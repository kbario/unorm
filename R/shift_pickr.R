#' Peak-Picking
#'
#' @param x
#' @param ppm
#' @param shift
#'
#' @return
#' @export
#'
#' @examples
shift_pickr <- function(X, ppm = NULL, search = c(3,3.1), pm = 0.005){
  if (is.null(dim(X)) & is.null(length(X))){
    stop("X is neither a spectrum or a matrix, please provide an approproiate X.")
  }
  if (is.null(dim(X)) & !is.null(length(X))){
    if (is.null(ppm)){
      stop('Please provide a X-matched ppm variable')
    }
    p <- ppm
    if (length(search)>2){
      i <- search
    } else {
      i <- get_idx(search, p)
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
    if (length(search)>2){
      i <- search
    } else {
      i <- get_idx(search, p)
    }
    s <- t(sapply(1:nrow(X), function(y){
      m <- unname(which(X[y,]==max(X[y,i])))
      s <- get_idx(c(p[m]-pm,p[m]+pm), p)
      return(s)
    }))
    return(s)
  }
}
